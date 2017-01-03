#include <deal.II/base/multithread_info.h>  
#include <deal.II/base/parameter_handler.h>

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>   
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream> 
#include <string>
#include <vector>
#include <DataTools.h>
namespace TRL
{
  using namespace dealii;
#include "InitialValue.h"
#include "parameters.h"

  template <int dim>
  class Heat_Pipe
  {
  public:
    Heat_Pipe(int argc, char *argv[]);
    ~Heat_Pipe();
    void run();

  private:
    void read_grid_temperature();
    void setup_system_temperature();
    void assemble_system_temperature();
    void solve_temperature();
    void initial_condition_temperature();

    void output_results ();
    void fill_output_vectors();
    void update_met_data ();

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;

    SparseMatrix<double> system_matrix;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix_new;
    SparseMatrix<double> laplace_matrix_old;
    Vector<double>       system_rhs;
    Vector<double>       solution;
    Vector<double>       old_solution;

    unsigned int timestep_number_max;
    unsigned int timestep_number;
    double       time;
    double       time_step;
    double       time_max;
    double       theta_temperature;

    Threads::Mutex assembler_lock;
    Parameters::AllParameters<dim>  parameters;

    std::vector< std::vector<int> >    date_and_time;
    std::vector< std::vector<double> > met_data;
    std::vector< std::vector<double> > depths_coordinates;
    std::vector< std::vector<double> > temperatures_at_points;
    std::vector< std::vector<double> > point_source_magnitudes;
    double old_room_temperature, new_room_temperature;
    double old_surface_temperature, new_surface_temperature;
    double old_point_source_magnitude, new_point_source_magnitude;

    std::ofstream output_file;
  };

  template<int dim>
  Heat_Pipe<dim>::Heat_Pipe(int argc, char *argv[])
    :
    dof_handler(triangulation),
    fe(1)
  {
    if (argc!=2)
      {
	std::cout << "Wrong number of input arguments.\n"
		  << "Number of arguments passed: " << argc << "\n"
		  << "Number of arguments expected: 2\n"
		  << "Missing input file?\n" << std::endl;
	throw 1;
      }

    std::string input_filename = argv[1];
    std::cout << "parameter file: " << input_filename << "\n";

    ParameterHandler prm;
    Parameters::AllParameters<dim>::declare_parameters (prm);
    prm.read_input(input_filename);
    parameters.parse_parameters (prm);

    theta_temperature   = parameters.theta;
    timestep_number_max = parameters.timestep_number_max;
    time_step           = parameters.time_step;
    time_max            = time_step*timestep_number_max;

    /*
      We want to read the file containing the coordinates
      of the points we are interested in. We will store
      them in a vector and use them later to extract data
      from the solution vector and do some calculations
      (e.g. stored thermal energy).
    */
    std::vector< std::vector<int> >    dummy_matrix;
    std::vector< std::string >         filenames;
    filenames.push_back(parameters.depths_file);

    DataTools data_tools;
    data_tools.read_data (filenames,
			  dummy_matrix,
			  depths_coordinates,
			  false);

    std::cout << "Available depth coordinate entries: "
	      << depths_coordinates.size() << std::endl
	      << "Depth coordinates (m):\n"
	      << "\tX\tY\tZ\n";
    for (unsigned int i=0; i<depths_coordinates.size(); i++)
      {
	for (unsigned int j=0; j<depths_coordinates[i].size(); j++)
	  std::cout << "\t" << depths_coordinates[i][j];
	std::cout << "\n";
      }

    std::string output_filename="output_data.txt";
    remove(output_filename.c_str());
    output_file.open(output_filename.c_str(),std::ios::app);
    if (!output_file.is_open()) //some error with the file
      {
	std::cout << "Error opening \"output_data.txt\" file\n";
	throw 1;
      }
    old_room_temperature       = 0.;
    new_room_temperature       = 0.;
    old_surface_temperature    = 0.;
    new_surface_temperature    = 0.;
    old_point_source_magnitude = 0.;
    new_point_source_magnitude = 0.;
    time=0.;
    timestep_number=0;
  }

  template<int dim>
  Heat_Pipe<dim>::~Heat_Pipe ()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void Heat_Pipe<dim>::read_grid_temperature()
  {
    GridGenerator::hyper_cube (triangulation,-1.*parameters.domain_size, 0);
    triangulation.refine_global (parameters.refinement_level);
    dof_handler.distribute_dofs (fe);
  }

  template <int dim>
  void Heat_Pipe<dim>::setup_system_temperature()
  {
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
					     hanging_node_constraints);
    hanging_node_constraints.close ();

    DynamicSparsityPattern csp(dof_handler.n_dofs(),
			       dof_handler.n_dofs());

    DoFTools::make_sparsity_pattern (dof_handler, csp);

    hanging_node_constraints.condense (csp);
    sparsity_pattern.copy_from (csp);
  }

  template <int dim>
  void Heat_Pipe<dim>::assemble_system_temperature()
  {
    system_rhs.reinit         (dof_handler.n_dofs());
    system_matrix.reinit      (sparsity_pattern);
    mass_matrix.reinit        (sparsity_pattern);
    laplace_matrix_new.reinit (sparsity_pattern);
    laplace_matrix_old.reinit (sparsity_pattern);

    //---------------------------------------------
    QGauss<dim>   quadrature_formula(3);
    const QGauss<dim-1>   face_quadrature_formula(3);
    FEValues<dim> fe_values(fe, quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
				     update_values | update_gradients |
				     update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
 
    FullMatrix<double> cell_mass_matrix        (dofs_per_cell,dofs_per_cell);
    FullMatrix<double> cell_laplace_matrix_new (dofs_per_cell,dofs_per_cell);
    FullMatrix<double> cell_laplace_matrix_old (dofs_per_cell,dofs_per_cell);
    Vector<double>     cell_rhs                (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	fe_values.reinit (cell);
	cell_mass_matrix        = 0;
	cell_laplace_matrix_new = 0;
	cell_laplace_matrix_old = 0;
	cell_rhs                = 0;

	/*
	 * We are assuming that each layer is composed of three fractions: solid, liquid, frozen liquid,
	 * and gas. At the moment, the assumption is that liquid is water, frozen liquid is ice and gas
	 * is air for all layers. Solids can vary, in out particular case, at least one layer has a
	 * different kind of solids.
	 *
	 * Porosity and degree of saturation is also layer dependent. Thus, fractions are layer
	 * dependent.
	 *
	 * Thermal properties can be considered equal for liquids and gas fractions but not for solids
	 * across the layers.
	 *
	 * NOTE that at the moment there are two layers composed of a single element (plastic lids)
	 * the thermal properties of these layers are calculated in a different way until we think of
	 * a way of homogenize the code.
	 * */

	double thermal_conductivity                 = -1.E10;
	double total_volumetric_heat_capacity       = -1.E10;
	double old_heat_loss                           = 0.;
	double new_heat_loss                           = 0.;
	{
	  double specific_heat_capacity_liquids = parameters.specific_heat_capacity_liquids;
	  double specific_heat_capacity_ice     = parameters.specific_heat_capacity_ice;
	  double specific_heat_capacity_gas     = parameters.specific_heat_capacity_air;
	  double density_liquids                = parameters.density_liquids;
	  double density_ice                    = parameters.density_ice;
	  double density_gas                    = parameters.density_air;
	  
	  double freezing_point                 = parameters.freezing_point;
	  double reference_temperature          = parameters.reference_temperature;
	  double coefficient_alpha              = parameters.alpha;
	  double latent_heat_of_fusion          = parameters.latent_heat;
	  
	  double cell_center=
	    cell->center()[0];
	  double cell_temperature=
	    0.5*VectorTools::point_value(dof_handler,solution    ,Point<dim>(cell->center()[0]))+
	    0.5*VectorTools::point_value(dof_handler,old_solution,Point<dim>(cell->center()[0]));
	  double convective_coefficient = 11.6516607672; // W/m3K
	  old_heat_loss = -1.*convective_coefficient*(cell_temperature-old_room_temperature); //W/m3
	  new_heat_loss = -1.*convective_coefficient*(cell_temperature-new_room_temperature); //W/m3
			
	  double degree_of_saturation_ice           =0.;
	  double derivative_degree_of_saturation_ice=0.;
	  if (cell_temperature<=freezing_point)
	    {
		  degree_of_saturation_ice=
				  1.-pow(1.-(cell_temperature-freezing_point),coefficient_alpha);
		  derivative_degree_of_saturation_ice=
				  coefficient_alpha*pow(1.-(cell_temperature-freezing_point),coefficient_alpha-1.);
	    }
	  
	  double specific_heat_capacity_solids=-1.E10;
	  double density_solids               =-1.E-10;
	  double porosity                     =-1.E10;
	  double degree_of_saturation         =-1.-10;
	  if (cell_center>-1.*(parameters.material_0_depth+parameters.material_0_thickness))/*layer 0*/
	    {
	      thermal_conductivity          = parameters.material_0_thermal_conductivity;
	      specific_heat_capacity_solids = parameters.material_0_specific_heat_capacity;
	      density_solids                = parameters.material_0_density;
	      porosity                      = parameters.material_0_porosity;
	      degree_of_saturation          = parameters.material_0_degree_of_saturation;
	    }
	  else if (cell_center<=-1.*parameters.material_1_depth &&
		   cell_center>-1.*(parameters.material_1_depth+parameters.material_1_thickness))/*layer 1*/
	    {
	      thermal_conductivity          = parameters.material_1_thermal_conductivity;
	      specific_heat_capacity_solids = parameters.material_1_specific_heat_capacity;
	      density_solids                = parameters.material_1_density;
	      porosity                      = parameters.material_1_porosity;
	      degree_of_saturation          = parameters.material_1_degree_of_saturation;
	    }
	  else if (cell_center<=-1.* parameters.material_2_depth &&
		   cell_center>-1.*(parameters.material_2_depth+parameters.material_2_thickness))/*layer 2*/
	    {
	      thermal_conductivity          = parameters.material_2_thermal_conductivity;
	      specific_heat_capacity_solids = parameters.material_2_specific_heat_capacity;
	      density_solids                = parameters.material_2_density;
	      porosity                      = parameters.material_2_porosity;
	      degree_of_saturation          = parameters.material_2_degree_of_saturation;
	    }
	  else if (cell_center<=-1.* parameters.material_3_depth &&
		   cell_center>-1.*(parameters.material_3_depth+parameters.material_3_thickness))/*layer 3*/
	    {
	      thermal_conductivity          = parameters.material_3_thermal_conductivity;
	      specific_heat_capacity_solids = parameters.material_3_specific_heat_capacity;
	      density_solids                = parameters.material_3_density;
	      porosity                      = parameters.material_3_porosity;
	      degree_of_saturation          = parameters.material_3_degree_of_saturation;
	    }
	  else if (cell_center<=-1.*parameters.material_4_depth)/*layer 4*/
	    {
	      thermal_conductivity          = parameters.material_4_thermal_conductivity;
	      specific_heat_capacity_solids = parameters.material_4_specific_heat_capacity;
	      density_solids                = parameters.material_4_density;
	      porosity                      = parameters.material_4_porosity;
	      degree_of_saturation          = parameters.material_4_degree_of_saturation;
	    }
	  else
	    {
	      std::cout << "Error. Cell center not found." << std::endl;
	      throw -1;
	    }

	  double Hc=
	    (1.-degree_of_saturation_ice)*
		porosity*degree_of_saturation*specific_heat_capacity_liquids*density_liquids
		+porosity*specific_heat_capacity_gas*density_gas*(1.-degree_of_saturation)
	    +specific_heat_capacity_solids*density_solids*(1.-porosity)
		+porosity*degree_of_saturation*degree_of_saturation_ice*specific_heat_capacity_ice*density_ice;

	  double a=
	    (cell_temperature-reference_temperature)*
	    (degree_of_saturation*density_ice*specific_heat_capacity_ice
	     -degree_of_saturation*density_liquids*specific_heat_capacity_liquids);
	  double b=
	    degree_of_saturation*density_ice*latent_heat_of_fusion;
	  
	  total_volumetric_heat_capacity=
	    (1.-degree_of_saturation_ice)*
		porosity*degree_of_saturation*specific_heat_capacity_liquids*density_liquids
		+porosity*specific_heat_capacity_gas*density_gas*(1.-degree_of_saturation)
		+specific_heat_capacity_solids*density_solids*(1.-porosity)
		+porosity*degree_of_saturation*degree_of_saturation_ice*specific_heat_capacity_ice*density_ice
	    +
		porosity*derivative_degree_of_saturation_ice*
	    ((cell_temperature-reference_temperature)*
	     (degree_of_saturation*density_ice*specific_heat_capacity_ice
	      -degree_of_saturation*density_liquids*specific_heat_capacity_liquids)
	     -degree_of_saturation*density_ice*latent_heat_of_fusion);

	  if (thermal_conductivity<0. || total_volumetric_heat_capacity<0.)
	    {
	      std::cout << "thermal_conductivity: " << thermal_conductivity << "\t"
			<< "total_volumetric_heat_capacity: "<< total_volumetric_heat_capacity << "\t"
			<< "cell temperature: " << cell_temperature << "\t"
			<< "ice content: " << degree_of_saturation_ice << "\t"
			<< "dSi/dT: " << derivative_degree_of_saturation_ice << "\t"
			<< "Hc: " << Hc << "\t" 
			<< "a: " << a << "\t"
			<< "b: " << b << "\t"
			<< "a-b: " << a-b << "\t" << std::endl;
	      throw -1;
	    }
	}

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  {
	    /*
	     * Here is were we assemble the matrices and vectors that appear after
	     * we discretize the problem in space and time using the finite element
	     * method. And here is also were we need to put any sinks or sources we
	     * want to implement. For the moment the magnitud of the source is user
	     * defined (by an external file). But at some point this must be
	     * calculated in a more appropiated way.
	     */
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    cell_mass_matrix(i,j)+=
		      total_volumetric_heat_capacity*
		      fe_values.shape_value(i,q_point) *
		      fe_values.shape_value(j,q_point) *
		      fe_values.JxW(q_point);
		    cell_laplace_matrix_new(i,j)+=(thermal_conductivity *
						   fe_values.shape_grad(i,q_point) *
						   fe_values.shape_grad(j,q_point) *
						   fe_values.JxW(q_point));
		    cell_laplace_matrix_old(i,j)+=(thermal_conductivity *
						   fe_values.shape_grad(i,q_point) *
						   fe_values.shape_grad(j,q_point) *
						   fe_values.JxW(q_point));
		  }
		cell_rhs(i)+=
		  new_heat_loss*theta_temperature*time_step*
		  fe_values.shape_value(i,q_point) *
		  fe_values.JxW(q_point)
		  +
		  old_heat_loss*(1-theta_temperature)*time_step*
		  fe_values.shape_value(i,q_point) *
		  fe_values.JxW(q_point);
	      }
	  }
	
	if (parameters.fixed_at_top==false)
	{
		for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			if (cell->face(face)->at_boundary() &&
					fabs(cell->face(face)->center()[0]+0.)<0.0001)
			{
				fe_face_values.reinit (cell, face);
				for (unsigned int q_face_point=0; q_face_point<n_face_q_points; ++q_face_point)
					for (unsigned int i=0; i<dofs_per_cell; ++i)
					{
						cell_rhs(i) -=
								100*time_step*theta_temperature*
								(fe_face_values.shape_value(i,q_face_point) *
										fe_face_values.JxW(q_face_point))
										+
										100*time_step*(1.-theta_temperature)*
										(fe_face_values.shape_value(i,q_face_point) *
												fe_face_values.JxW(q_face_point));
					}
			}
	}
	cell->get_dof_indices (local_dof_indices);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		  laplace_matrix_new.add (local_dof_indices[i],local_dof_indices[j],cell_laplace_matrix_new(i,j));
			laplace_matrix_old.add (local_dof_indices[i],local_dof_indices[j],cell_laplace_matrix_old(i,j));
			mass_matrix.add        (local_dof_indices[i],local_dof_indices[j],cell_mass_matrix(i,j)       );
		}
		system_rhs(local_dof_indices[i]) += cell_rhs(i);
	}
      }
    
    
    Vector<double> tmp      (solution.size ());
    
    /*
      This is the section where the point source is included.
      Notice that there are other ways to do this, but the
      library has a function that is intended exactly for this
      kind of problem. Check the documentation of
      'create_point_source_vector' in deal.ii
    */
    if (parameters.point_source==true)
      {
	Point<dim> p(-1.*parameters.point_source_depth);
	VectorTools::create_point_source_vector(dof_handler,
						p,
						tmp);
	system_rhs.add           (old_point_source_magnitude* // (W/m3)
				  (1-theta_temperature)*time_step
				  +
				  new_point_source_magnitude* // (W/m3)
				  (  theta_temperature)*time_step,tmp);
      }
    //--------------------------------------------------------
    mass_matrix.vmult        ( tmp,old_solution);
    system_rhs.add           ( 1.0,tmp);
    laplace_matrix_old.vmult ( tmp,old_solution);
    system_rhs.add           (-(1 - theta_temperature) * time_step,tmp);
    
    system_matrix.copy_from (mass_matrix);
    system_matrix.add       (theta_temperature * time_step, laplace_matrix_new);
    
    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);
    
    if (parameters.fixed_at_bottom)
      {
	std::map<unsigned int,double> boundary_values;
	
	VectorTools::interpolate_boundary_values (dof_handler,
						  0,
						  ConstantFunction<dim>(parameters.bottom_fixed_value),
						  boundary_values);
	MatrixTools::apply_boundary_values (boundary_values,
					    system_matrix,
					    solution,
					    system_rhs);
      }
    if (parameters.fixed_at_top)
      {
	std::map<unsigned int,double> boundary_values;
	VectorTools::interpolate_boundary_values (dof_handler,
						  1,
						  ConstantFunction<dim>(parameters.theta * new_surface_temperature +
									(1-parameters.theta) * old_surface_temperature),
						  boundary_values);
	MatrixTools::apply_boundary_values (boundary_values,
					    system_matrix,
					    solution,
					    system_rhs);
      }
  }

  template <int dim>
  void Heat_Pipe<dim>::solve_temperature()
  {
    SolverControl solver_control (solution.size(),
				  1e-8*system_rhs.l2_norm ());
    SolverCG<> cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize (system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs,
	      preconditioner);

    hanging_node_constraints.distribute (solution);
  }

  template <int dim>
  void Heat_Pipe<dim>::output_results()
  {
    
	  DataOut<dim> data_out;

	  data_out.attach_dof_handler(dof_handler);
	  data_out.add_data_vector(solution,"solution");
	  data_out.build_patches();

	  std::stringstream t;
	  t << timestep_number;

	  std::stringstream d;
	  d << dim;

	  std::string filename = parameters.output_directory + "/solution_"
			  + d.str() + "d_time_"
			  + t.str() + ".vtu";

	  std::ofstream output (filename.c_str());
	  data_out.write_vtu (output);

	  // extract and save temperatures
	  // at selected coordinates
	  if (dim==1)
	  {
	    /*
	      Extract temperatures from the solution vector
	    */
		  std::vector<double> temp_vector;
		  if (temperatures_at_points.size()==0)
		  {
			  for (unsigned int i=0; i<depths_coordinates.size(); i++)
			    temp_vector
			      .push_back(VectorTools::point_value(dof_handler,solution,
								  Point<dim>(-1.*depths_coordinates[i][2])));

		  }
		  else
		  {
			  for (unsigned int i=0; i<depths_coordinates.size(); i++)
				  temp_vector
				  .push_back(VectorTools::point_value(dof_handler,solution,
						  Point<dim>(-1.*depths_coordinates[i][2])));
		  }
		  temperatures_at_points.push_back(temp_vector);
		  /*
	  Save them to some file.
		   */
		  output_file << timestep_number;
		  for (unsigned int i=0; i<temp_vector.size(); i++)
		    output_file << "\t" << std::setprecision(5) << temp_vector[i];
		  //output_file << std::endl;


		  //	Estimate thermal energy

		  double thermal_energy=0.;
		  typename DoFHandler<dim>::active_cell_iterator
		  cell = dof_handler.begin_active(),
		  endc = dof_handler.end();
		  for (; cell!=endc; ++cell)
		  {
			  double volumetric_heat_capacity       = 0.;
			  double cell_center                    = cell->center()[0];
			  double cell_temperature               =
					  0.5*VectorTools::point_value(dof_handler,solution    ,Point<dim>(cell->center()[0]))+
					  0.5*VectorTools::point_value(dof_handler,old_solution,Point<dim>(cell->center()[0]));

			  double specific_heat_capacity_solids  = -1.E-10;
			  double density_solids                 = -1.E-10;
			  double porosity                       = -1.E-10;
			  double degree_of_saturation           = -1.E-10;
			  double specific_heat_capacity_liquids = parameters.specific_heat_capacity_liquids;
			  double specific_heat_capacity_ice     = parameters.specific_heat_capacity_ice;
			  double specific_heat_capacity_gas     = parameters.specific_heat_capacity_air;
			  double density_liquids                = parameters.density_liquids;
			  double density_ice                    = parameters.density_ice;
			  double density_gas                    = parameters.density_air;
			  double freezing_point                 = parameters.freezing_point;
			  double reference_temperature          = parameters.reference_temperature;
			  double coefficient_alpha              = parameters.alpha;
			  double latent_heat_of_fusion          = parameters.latent_heat;
			  double degree_of_saturation_ice       = 0. ;
			  if (cell_temperature<=freezing_point)
			  {
				  degree_of_saturation_ice=
						  1.-pow(1.-(cell_temperature-freezing_point),coefficient_alpha);
			  }
			  if (cell_center>-1.*(parameters.material_0_depth+parameters.material_0_thickness))       /*layer 0*/
			  {
				  specific_heat_capacity_solids = parameters.material_0_specific_heat_capacity;
				  density_solids                = parameters.material_0_density;
				  porosity                      = parameters.material_0_porosity;
				  degree_of_saturation          = parameters.material_0_degree_of_saturation;
			  }
			  else if (cell_center<=-1.*parameters.material_1_depth &&
					  cell_center>-1.*(parameters.material_1_depth+parameters.material_1_thickness))  /*layer 1*/
			  {
				  specific_heat_capacity_solids = parameters.material_1_specific_heat_capacity;
				  density_solids                = parameters.material_1_density;
				  porosity                      = parameters.material_1_porosity;
				  degree_of_saturation          = parameters.material_1_degree_of_saturation;
			  }
			  else if (cell_center<=-1.* parameters.material_2_depth &&
					  cell_center>-1.*(parameters.material_2_depth+parameters.material_2_thickness))  /*layer 2*/
			  {
				  specific_heat_capacity_solids = parameters.material_2_specific_heat_capacity;
				  density_solids                = parameters.material_2_density;
				  porosity                      = parameters.material_2_porosity;
				  degree_of_saturation          = parameters.material_2_degree_of_saturation;
			  }
			  else if (cell_center<=-1.* parameters.material_3_depth &&
					  cell_center>-1.*(parameters.material_3_depth+parameters.material_3_thickness))  /*layer 3*/
			  {
				  specific_heat_capacity_solids = parameters.material_3_specific_heat_capacity;
				  density_solids                = parameters.material_3_density;
				  porosity                      = parameters.material_3_porosity;
				  degree_of_saturation          = parameters.material_3_degree_of_saturation;
			  }
			  else if (cell_center<=-1.*parameters.material_4_depth)                                   /*layer 4*/
			  {
				  specific_heat_capacity_solids = parameters.material_4_specific_heat_capacity;
				  density_solids                = parameters.material_4_density;
				  porosity                      = parameters.material_4_porosity;
				  degree_of_saturation          = parameters.material_4_degree_of_saturation;
			  }
			  else
			  {
				  std::cout << "Error. Thermal energy wrong calculation." << std::endl;
				  throw -1;
			  }
			  volumetric_heat_capacity=
					  (1.-degree_of_saturation_ice)*
					  porosity*degree_of_saturation*specific_heat_capacity_liquids*density_liquids
					  +porosity*specific_heat_capacity_gas*density_gas*(1.-degree_of_saturation)
					  +specific_heat_capacity_solids*density_solids*(1.-porosity)
					  +porosity*degree_of_saturation*degree_of_saturation_ice*specific_heat_capacity_ice*density_ice;

			  thermal_energy+=
					  cell->diameter()*
					  (volumetric_heat_capacity*(/*cell_temperature*/VectorTools::point_value(dof_handler,solution    ,Point<dim>(cell->center()[0]))-reference_temperature)
							  -latent_heat_of_fusion*porosity*degree_of_saturation*degree_of_saturation_ice*density_ice);
		  }
		  /*
	  Add the estimated energy to the output file
		   */
		  output_file << "\t" << std::setprecision(5) << thermal_energy;
		  output_file << std::endl;
		  //output_file.close();
		  std::cout << "Estimated thermal energy: " << thermal_energy << "\n";
	  }
	  else
	  {
		  std::cout << "Error in output_results function.\n"
				  << "Currently implement only for 1D\n";
		  throw 1;
	  }
  }

  template <int dim>
  void Heat_Pipe<dim>::update_met_data ()
  {
    /*
      Originally this function read a file with date (dd/mm/yyyy) and met data (air
      temperature, solar radiation, wind speed, etc) (hence the names). For the
      purpose of analysing a simplified composite region in 1D with fixed top boundary
      condition I'm assuming that we are reading a file with a single column corresponding
      to surface temperature at every time step. We can come back to more complex met
      data files later on
    */
    if (met_data.size()==0)
      {
	DataTools data_tools;
	std::vector<std::string> filenames;
	filenames.push_back(parameters.top_fixed_value_file);
	data_tools.read_data (filenames,
			      date_and_time,
			      met_data,
			      false);

	std::cout << "\tAvailable surface data lines: " << met_data.size()
		  << std::endl << std::endl;

	if (parameters.point_source==true)
	  {
	    std::vector< std::vector<int> > dummy;
	    filenames.clear();
	    filenames.push_back(parameters.point_source_file);
	    data_tools.read_data (filenames,
				  dummy,
				  point_source_magnitudes,
				  false);

	    std::cout << "\tAvailable point source entries: "
		      << point_source_magnitudes.size()
		      << std::endl;
	  }
      }
    old_room_temperature    = met_data[timestep_number-1][1];
    new_room_temperature    = met_data[timestep_number  ][1];
    old_surface_temperature = met_data[timestep_number-1][0];
    new_surface_temperature = met_data[timestep_number  ][0];

    if (parameters.point_source==true)
      {
	old_point_source_magnitude = point_source_magnitudes [timestep_number-1][1];
	new_point_source_magnitude = point_source_magnitudes [timestep_number  ][1];
      }
  }

  template <int dim>
  void Heat_Pipe<dim>::initial_condition_temperature()
  {
    /*
      Here the vectors containing the name of the file with 
      the initial condition is defined
    */
    std::vector< std::string > filenames;
    filenames.push_back(parameters.initial_condition_file);
    /*
      And the vector containing the actual data (depth and
      temperature). Note that it is actually a vector of
      vector so, rather, a matrix.
      We use an external function to read the data and fill
      the initial condition matrix.
    */
    std::vector< std::vector<int> > dummy_matrix;
    std::vector< std::vector<double> > initial_condition;
    DataTools data_tools;
    data_tools.read_data (filenames,
			  dummy_matrix,
			  initial_condition,
			  false);
    // we use the following lines to print out the vectors 
    // for (unsigned int i=0; i<initial_condition.size(); ++i) 
    //   for (unsigned int j=0; j<initial_condition[i].size(); ++j)
    // 	std::cout << initial_condition[i][j] << "\n";
    

    /*
      The format of the matrix need to be changed from:

      ' std::vector< std::vector<double> > '
      to
      ' std::vector< std::pair<double,double> > '

      this is because the function that interpolates the 
      depths from the provided file takes this kind of
      format. Note that it is assumed that the matrix
      has two elements per row. If more are provided
      they will be ignored.
    */
    std::vector< std::pair<double,double> > initial_condition_table;
    for (unsigned int i=0; i<initial_condition.size(); i++)
      initial_condition_table
	.push_back(std::make_pair(initial_condition[i][0],
				  initial_condition[i][1]));
    /*
      Print number of lines available in the initial
      condition file and the actual data read.
      Currently it is expected a file with the following 
      format:
      Depth (m) \t Temperature (C)  <-- this line is not expected in the file

      0.0             T0
      d1              T1
      ...             ...
      dN              TN

      Note that the file starts with the temperature at x=0 m.
      If other depth is provided in the first line, the program
      will probalbly throw a not very informative error. Also,
      all depth values are positive.
    */
    std::cout << "Available initial condition entries: "
	      << initial_condition.size()    << std::endl
	      << "Initial condition: \n\tDepth\tTemperature (C)\n";
    for (unsigned int i=0; i<initial_condition.size(); i++)
      std::cout << "\t" << initial_condition[i][0] << "\t" << initial_condition[i][1] <<std::endl; 

    VectorTools::project (dof_handler,
			  hanging_node_constraints,
			  QGauss<dim>(2),
			  InitialValue<dim>(initial_condition_table),
			  old_solution);

    // VectorTools::project (dof_handler,
    // 			  hanging_node_constraints,
    // 			  QGauss<dim>(3),
    // 			  InitialValue<dim>(10.),
    // 			  old_solution);
    solution=old_solution;
  }


  
  template <int dim>
  void Heat_Pipe<dim>::run()
  {
	  read_grid_temperature();
	  setup_system_temperature();
	  solution.reinit (dof_handler.n_dofs());
	  old_solution.reinit (dof_handler.n_dofs());
	  initial_condition_temperature();

	  {
		  std::cout.setf( std::ios::fixed);
		  std::cout.precision(3);
		  std::cout << "\tPosition of material layers:\n"
				  << "\t\tLayer 1: from "
				  << parameters.material_0_depth << " to "
				  << parameters.material_0_depth+parameters.material_0_thickness << "\t"
				  << "k:" << parameters.material_0_thermal_conductivity << " W/mK\t"
				  << "Cp:" << (parameters.material_0_specific_heat_capacity*
						  parameters.material_0_density)/1000. << " MJ/m3K\n"
						  << "\t\tLayer 2: from "
						  << parameters.material_1_depth << " to "
						  << parameters.material_1_depth+parameters.material_1_thickness << "\t"
						  << "k:" << parameters.material_1_thermal_conductivity << " W/mK\t"
						  << "Cp:" << (parameters.material_1_specific_heat_capacity*
								  parameters.material_1_density)/1000. << " MJ/m3K\n"
								  << "\t\tLayer 3: from "
								  << parameters.material_2_depth << " to "
								  << parameters.material_2_depth+parameters.material_2_thickness << "\t"
								  << "k:" << parameters.material_2_thermal_conductivity << " W/mK\t"
								  << "Cp:" << (parameters.material_2_specific_heat_capacity*
										  parameters.material_2_density)/1000. << "MJ/m3K\n"
										  << "\t\tLayer 4: from "
										  << parameters.material_3_depth << " to "
										  << parameters.material_3_depth+parameters.material_3_thickness << "\t"
										  << "k:" << parameters.material_3_thermal_conductivity << " W/mK\t"
										  << "Cp:" << (parameters.material_3_specific_heat_capacity*
												  parameters.material_3_density)/1000. << "MJ/m3K\n"
												  << "\t\tLayer 5: from "
												  << parameters.material_4_depth << " to "
												  << parameters.material_4_depth+parameters.material_4_thickness << "\t"
												  << "k:" << parameters.material_4_thermal_conductivity << " W/mK\t"
												  << "Cp:" << (parameters.material_4_specific_heat_capacity*
														  parameters.material_4_density)/1000. << "MJ/m3K\n";
	  }

	  for (time=time_step, timestep_number=1;
			  time<=time_max;
			  time+=time_step, ++timestep_number)
	  {
	    if (parameters.fixed_at_top) // To correct the error in calculation when we change the B.C. between 1st type and 2nd type problem
	      update_met_data();
	    
		  int iteration=0;
		  double total_error =1E10;
		  double solution_l1_norm_previous_iteration;
		  double solution_l1_norm_current_iteration;
		  do
		  {
			  assemble_system_temperature();

			  solution_l1_norm_previous_iteration=solution.l2_norm();
			  solve_temperature();
			  solution_l1_norm_current_iteration=solution.l2_norm();
			  
			  total_error=
			    1.-std::fabs(solution_l1_norm_previous_iteration/solution_l1_norm_current_iteration);
			  iteration++;
			  
			  std::cout << "\ttime step: " << timestep_number << "\tTo" << old_surface_temperature << "\tTn" << new_surface_temperature << "\n";
		  }while (std::fabs(total_error)>5E-3);
		  
		  std::cout << "Time step " << timestep_number << "\titeration: " << iteration << "\ttotal_error: " << total_error  << "\n";
		  if (parameters.output_frequency!=0 && timestep_number%parameters.output_frequency==0)
		    {
		      output_results();
		    }
		  
		  old_solution=solution;
	  }
	  output_file.close();
	  std::cout << "\t Job Done!!"
			  << std::endl;
  }
}

int main (int argc, char *argv[])
{
  try
    {
      using namespace TRL;
      using namespace dealii;
      {
	deallog.depth_console (0);

	Heat_Pipe<1> laplace_problem(argc,argv);

	laplace_problem.run();
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
