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

#include <PorousMaterial.h>
#include <DataTools.h>
#include <Names.h>

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream> 
#include <string>
#include <vector>
#include <tuple>

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

    void material_data(const double cell_center/*(m)*/,
		       const double cell_temperature/*(C)*/,
		       double &thermal_conductivity/*(W/mK)*/,
		       double &total_volumetric_heat_capacity/*(J/m3K)*/,
		       double &ice_saturation);
    double cell_thermal_energy(const double cell_center/*(m)*/,
			       const double cell_temperature/*(C)*/,
			       const double cell_diameter/*(m)*/);
    double thermal_losses(const double temperature_gradient/*(m)*/);
    unsigned int find_layer(double cell_center);
    //double snow_surface_heat_flux(double surface_temperature); //(W/m2)

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

    //std::vector< std::vector<int> >    date_and_time;
    std::vector< std::vector<double> > met_data;
    std::vector<std::vector<double> > interpolated_temperature_surface;
    std::vector<std::vector<double> > interpolated_temperature_room;
    std::vector< std::vector<double> > depths_coordinates;
    std::vector< std::vector<double> > temperatures_at_points;
    std::vector< std::vector<double> > point_source_magnitudes;
    double old_room_temperature, new_room_temperature;
    double old_surface_temperature, new_surface_temperature;
    double old_point_source_magnitude, new_point_source_magnitude;
    double column_thermal_energy;
    double thermal_conductivity_liquids;
    double thermal_conductivity_air;

    std::ofstream output_file;
    /*
     * string "material_name"
     * double "porosity"
     * double "degree_of_saturation"
     * string "relationship"
     */
    std::vector<std::tuple<std::string,double,double,std::string> > layer_data;
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

    const std::string input_filename = argv[1];
    std::cout << "parameter file: " << input_filename << "\n";

    ParameterHandler prm;
    Parameters::AllParameters<dim>::declare_parameters (prm);

    std::ifstream inFile;
    inFile.open(input_filename.c_str());
    
    //prm.read_input(input_filename);
    prm.parse_input(inFile,input_filename);
    
    parameters.parse_parameters (prm);

    theta_temperature   = parameters.theta;
    timestep_number_max = parameters.timestep_number_max;
    time_step           = parameters.time_step;
    time_max            = time_step*timestep_number_max;

    thermal_conductivity_liquids   = parameters.thermal_conductivity_liquids;
    thermal_conductivity_air       = parameters.thermal_conductivity_air;
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
			  depths_coordinates);

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

    std::string output_filename=parameters.output_file;
    remove(output_filename.c_str());
    output_file.open(output_filename.c_str(),std::ios::app);
    if (!output_file.is_open()) //some error with the file
      {
    	std::cout << "Error opening output data file\n";
    	throw 1;
      }

    if (parameters.boundary_condition_top.compare("first")!=0 &&
    	parameters.boundary_condition_top.compare("second")!=0 &&
	parameters.boundary_condition_top.compare("third")!=0 )
      {
    	std::cout << "\n\n\tError, wrong top boundary condition type\n\n";
      }

    old_room_temperature       = 0.;
    new_room_temperature       = 0.;
    old_surface_temperature    = 0.;
    new_surface_temperature    = 0.;
    old_point_source_magnitude = 0.;
    new_point_source_magnitude = 0.;
    time=0.;
    timestep_number=0;
    column_thermal_energy=0.;

    layer_data
      .push_back(std::tuple<std::string,double,double,std::string>
		 (parameters.material_0_name,
		  parameters.material_0_porosity,
		  parameters.material_0_degree_of_saturation,
		  parameters.material_0_thermal_conductivity_relationship));
    layer_data
      .push_back(std::tuple<std::string,double,double,std::string>
		 (parameters.material_1_name,
		  parameters.material_1_porosity,
		  parameters.material_1_degree_of_saturation,
		  parameters.material_1_thermal_conductivity_relationship));
    layer_data
      .push_back(std::tuple<std::string,double,double,std::string>
		 (parameters.material_2_name,
		  parameters.material_2_porosity,
		  parameters.material_2_degree_of_saturation,
		  parameters.material_2_thermal_conductivity_relationship));
    layer_data
      .push_back(std::tuple<std::string,double,double,std::string>
		 (parameters.material_3_name,
		  parameters.material_3_porosity,
		  parameters.material_3_degree_of_saturation,
		  parameters.material_3_thermal_conductivity_relationship));
     layer_data
      .push_back(std::tuple<std::string,double,double,std::string>
		 (parameters.material_4_name,
		  parameters.material_4_porosity,
		  parameters.material_4_degree_of_saturation,
		  parameters.material_4_thermal_conductivity_relationship));
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
  void Heat_Pipe<dim>::material_data(const double cell_center,
				     const double cell_temperature,
				     double &thermal_conductivity,
				     double &total_volumetric_heat_capacity,
				     double &ice_saturation)
  {
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
     *
     * The function receives thermal conductivity of solids in the layer of interest.
     * */
    /*
     * These variables are assumed to be constants. That's why we defined them inside the function.
     * */
    unsigned int layer_number=
      find_layer(cell_center);
    std::string material_name  =std::get<0>(layer_data[layer_number]);
    double porosity            =std::get<1>(layer_data[layer_number]);
    double degree_of_saturation=std::get<2>(layer_data[layer_number]);
    std::string relationship   =std::get<3>(layer_data[layer_number]);
    
    PorousMaterial porous_material(material_name,
				   porosity,
				   degree_of_saturation);
    thermal_conductivity=
      porous_material.thermal_conductivity(relationship);
    total_volumetric_heat_capacity=
      porous_material.volumetric_heat_capacity(cell_temperature);
    ice_saturation=
      porous_material.degree_of_saturation_ice(cell_temperature);
    
    if (thermal_conductivity<0. || total_volumetric_heat_capacity<0.)
      {
	std::cout << "thermal_conductivity: " << thermal_conductivity << "\t"
		  << "total_volumetric_heat_capacity: "<< total_volumetric_heat_capacity << "\t"
		  << "cell temperature: " << cell_temperature << "\t"
		  << "ice content: " << ice_saturation << "\t"
		  << std::endl;
	throw -1;
      }
  }

  template <int dim>
  double Heat_Pipe<dim>::cell_thermal_energy(const double cell_center,
					     const double cell_temperature,
					     const double cell_diameter)
  {
    unsigned int layer_number
      =find_layer(cell_center);
    std::string material_name  =std::get<0>(layer_data[layer_number]);
    double porosity            =std::get<1>(layer_data[layer_number]);
    double degree_of_saturation=std::get<2>(layer_data[layer_number]);
    
    PorousMaterial porous_material(material_name,
				   porosity,
				   degree_of_saturation);
    return cell_diameter*porous_material.thermal_energy(cell_temperature);
  }
  
  template <int dim>
  unsigned int Heat_Pipe<dim>::find_layer(double cell_center)
  {
    unsigned int layer_number=0;
    if (cell_center>-1.*(parameters.material_0_depth+parameters.material_0_thickness))
      layer_number=0;
    else if (cell_center<=-1.*parameters.material_1_depth &&
	     cell_center>-1.*(parameters.material_1_depth+parameters.material_1_thickness))
      layer_number=1;
    else if (cell_center<=-1.* parameters.material_2_depth &&
	     cell_center>-1.*(parameters.material_2_depth+parameters.material_2_thickness))
      layer_number=2;
    else if (cell_center<=-1.* parameters.material_3_depth &&
	     cell_center>-1.*(parameters.material_3_depth+parameters.material_3_thickness))
      layer_number=3;
    else if (cell_center<=-1.*parameters.material_4_depth)
      layer_number=4;
    else
      {
	std::cout << "Error. Cell centre not found." << std::endl;
	throw -1;
      }
    return layer_number;
  }
  
  template <int dim>
  double Heat_Pipe<dim>::thermal_losses(const double temperature_gradient)
  {
    /*
     * At the moment is the convective coefficient is read from the input
     * file, but this function could include the equations used to estimate
     * this coefficient.
     * */
    double convective_coefficient = parameters.heat_loss_factor; // W/m3K
    return (-1.*convective_coefficient*(temperature_gradient)); //(W/m3)
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

    QGauss<dim>   quadrature_formula(3);
    const QGauss<dim-1>   face_quadrature_formula(3);
    FEValues<dim> fe_values(fe, quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
				     update_values | update_gradients | update_normal_vectors|
				     update_quadrature_points | update_JxW_values);
	  
    const unsigned int dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
	  
    FullMatrix<double> cell_mass_matrix        (dofs_per_cell,dofs_per_cell);
    FullMatrix<double> cell_laplace_matrix_new (dofs_per_cell,dofs_per_cell);
    FullMatrix<double> cell_laplace_matrix_old (dofs_per_cell,dofs_per_cell);
    Vector<double>     cell_rhs                (dofs_per_cell);

    Vector<double> old_temperature_values (dofs_per_cell);
    Vector<double> new_temperature_values (dofs_per_cell);
	  
    std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
    std::vector<double> old_function_values     (n_q_points);
    std::vector<double> new_function_values     (n_q_points);
    column_thermal_energy= 0.;
    double flux_top=0.;
    double flux_bottom=0.;
    
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	cell_mass_matrix        = 0;
	cell_laplace_matrix_new = 0;
	cell_laplace_matrix_old = 0;
	cell_rhs                = 0;
	fe_values.reinit (cell);
	fe_values.get_function_values(old_solution,old_function_values);
	fe_values.get_function_values(    solution,new_function_values);

	new_temperature_values.reinit(cell->get_fe().dofs_per_cell);
	old_temperature_values.reinit(cell->get_fe().dofs_per_cell);

	cell->get_dof_values(solution    ,new_temperature_values);
	cell->get_dof_values(old_solution,old_temperature_values);
	
	double cell_thermal_conductivity          = -1.E10;
	double cell_total_volumetric_heat_capacity= -1.E10;
	double cell_ice_saturation                = -1.E10;
	
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  {
	    double average_cell_temperature=
	      (   theta_temperature)*new_function_values[q_point]+
	      (1.-theta_temperature)*old_function_values[q_point];
	    
	    double old_cell_heat_loss = thermal_losses(average_cell_temperature-old_room_temperature);
	    double new_cell_heat_loss = thermal_losses(average_cell_temperature-new_room_temperature);
	    
	    material_data(cell->center()[0],average_cell_temperature,
			  cell_thermal_conductivity,cell_total_volumetric_heat_capacity,
			  cell_ice_saturation);
	    
	    column_thermal_energy+=
	      cell_thermal_energy(cell->center()[0],average_cell_temperature,cell->diameter());
	    /*
	     * Here is were we assemble the matrices and vectors that appear after
	     * we discretize the problem in space and time using the finite element
	     * method. And here is also were we need to put any sinks or sources we
	     * want to implement. For the moment the magnitude of the source is user
	     * defined (by an external file). But at some point this must be
	     * calculated in a more appropiated way.
	     */
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    cell_mass_matrix(i,j)+=
		      cell_total_volumetric_heat_capacity*
		      fe_values.shape_value(i,q_point) *
		      fe_values.shape_value(j,q_point) *
		      fe_values.JxW(q_point);
		    cell_laplace_matrix_new(i,j)+=
		      cell_thermal_conductivity *
		      fe_values.shape_grad(i,q_point) *
		      fe_values.shape_grad(j,q_point) *
		      fe_values.JxW(q_point);
		    cell_laplace_matrix_old(i,j)+=
		      cell_thermal_conductivity *
		      fe_values.shape_grad(i,q_point) *
		      fe_values.shape_grad(j,q_point) *
		      fe_values.JxW(q_point);
		  }
		cell_rhs(i)+=
		  new_cell_heat_loss*theta_temperature*time_step*
		  fe_values.shape_value(i,q_point) *
		  fe_values.JxW(q_point)
		  +
		  old_cell_heat_loss*(1.-theta_temperature)*time_step*
		  fe_values.shape_value(i,q_point) *
		  fe_values.JxW(q_point);
	      }
	  }

	if (parameters.boundary_condition_top.compare("second")==0 ||
	    parameters.boundary_condition_top.compare("third")==0)
	  {
	    double convective_coefficient=10.;
	    double top_outbound_convective_coefficient=0.;
	    double top_inbound_heat_flux_new=0.;
	    double top_inbound_heat_flux_old=0.;
	    if (parameters.boundary_condition_top.compare("second")==0)
	      {
		top_inbound_heat_flux_new=-100.;
		top_inbound_heat_flux_old=-100.;
	      }
	    else if (parameters.boundary_condition_top.compare("third")==0)
	      {
		top_outbound_convective_coefficient=convective_coefficient;
		top_inbound_heat_flux_new=convective_coefficient*new_surface_temperature;
		top_inbound_heat_flux_old=convective_coefficient*old_surface_temperature;
	      }

	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	      if (cell->face(face)->at_boundary() &&
		  fabs(cell->face(face)->center()[0]+0.)<0.0001)
		{
		  fe_face_values.reinit (cell, face);
		  for (unsigned int q_face_point=0; q_face_point<n_face_q_points; ++q_face_point)
		    {
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  if(parameters.boundary_condition_top.compare("third")==0)
			    {
			      for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
				  cell_laplace_matrix_new (i,j)+=
				    top_outbound_convective_coefficient *
				    fe_face_values.shape_value (i,q_face_point) *
				    fe_face_values.shape_value (j,q_face_point) *
				    fe_face_values.JxW         (q_face_point);
				  cell_laplace_matrix_old (i,j)+=
				    top_outbound_convective_coefficient *
				    fe_face_values.shape_value (i,q_face_point) *
				    fe_face_values.shape_value (j,q_face_point) *
				    fe_face_values.JxW         (q_face_point);
				}
			    }
			  cell_rhs(i)+=
			    top_inbound_heat_flux_new*
			    time_step*theta_temperature*
			    fe_face_values.shape_value(i,q_face_point) *
			    fe_face_values.JxW(q_face_point)
			    +
			    top_inbound_heat_flux_old*
			    time_step*(1.-theta_temperature)*
			    fe_face_values.shape_value(i,q_face_point) *
			    fe_face_values.JxW(q_face_point);
			  
			  flux_top+=
			    top_inbound_heat_flux_new*
			    theta_temperature*
			    fe_face_values.shape_value(i,q_face_point) *
			    fe_face_values.JxW(q_face_point)
			    +
			    top_inbound_heat_flux_old*
			    (1.-theta_temperature)*
			    fe_face_values.shape_value(i,q_face_point) *
			    fe_face_values.JxW(q_face_point);
			}
		    }
		}
	      else if (cell->face(face)->at_boundary() &&
		       fabs(cell->face(face)->center()[0]+parameters.domain_size)<0.0001)
		{
		  fe_face_values.reinit (cell, face);
		  for (unsigned int q_face_point=0; q_face_point<n_face_q_points; ++q_face_point)
		    {
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    {			  
			      flux_bottom+=
				theta_temperature*
				cell_thermal_conductivity*
				fe_face_values.normal_vector(q_face_point)*
				fe_face_values.shape_grad (i,q_face_point)*
				new_temperature_values(i)*
				fe_face_values.shape_value(j,q_face_point)*
				fe_face_values.JxW(q_face_point)
				+
				(1.-theta_temperature)*
				cell_thermal_conductivity*
				fe_face_values.normal_vector(q_face_point)*
				fe_face_values.shape_grad (i,q_face_point)*
				old_temperature_values(i)*
				fe_face_values.shape_value(j,q_face_point) *
				fe_face_values.JxW(q_face_point);			      
			    }
			}
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

    std::cout << "\tflux top: " << flux_top << "\tflux bottom: " << flux_bottom << "\n";
    Vector<double> tmp(solution.size ());
    /*
     * This is the section where the point source is included.
     * Notice that there are other ways to do this, but the
     * library has a function that is intended exactly for this
     * kind of problem. Check the documentation of
     * 'create_point_source_vector' in deal.ii
     * */
    if (parameters.point_source==true)
      {
	Point<dim> p(-1.*parameters.point_source_depth);
	VectorTools::create_point_source_vector(dof_handler,p,tmp);
	system_rhs.add(old_point_source_magnitude*(1.-theta_temperature)*time_step
		       +new_point_source_magnitude*(   theta_temperature)*time_step,tmp);// (W/m3)
      }
    //--------------------------------------------------------
    mass_matrix.vmult        ( tmp,old_solution);
    system_rhs.add           ( 1.0,tmp);
    laplace_matrix_old.vmult ( tmp,old_solution);
    system_rhs.add           (-(1.-theta_temperature) * time_step,tmp);

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
    if (parameters.boundary_condition_top.compare("first")==0)
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
  void Heat_Pipe<dim>::fill_output_vectors()
  {
    /*
     * Extract and save temperatures at selected coordinates.
     **/
    if (dim==1)
      {
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
	//temp_vector.push_back(solution.l1_norm());
	temperatures_at_points.push_back(temp_vector);
	/*
	 * Save them to some file.
	 */
	output_file << timestep_number << "\t" << timestep_number*time_step;
	for (unsigned int i=0; i<temp_vector.size(); i++)
	  output_file << "\t" << std::setprecision(5) << temp_vector[i];

	output_file << "\t" << std::setprecision(5) << column_thermal_energy;
	output_file << std::endl;
      }
    else
      {
	std::cout << "Error in output_results function.\n"
		  << "Currently implement only for 1D\n";
	throw 1;
      }
  }

  template <int dim>
  void Heat_Pipe<dim>::output_results()
  {
    QGauss<dim>   quadrature_formula(3);
    FEValues<dim> fe_values(fe, quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);
    const unsigned int n_q_points      = quadrature_formula.size();
    std::vector<double> old_function_values     (n_q_points);
    std::vector<double> new_function_values     (n_q_points);    
    
    std::vector<double> ice_saturation_int;
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	fe_values.reinit (cell);
	fe_values.get_function_values(old_solution,old_function_values);
	fe_values.get_function_values(    solution,new_function_values);
	
	double average_cell_temperature=0.;
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  {
	    average_cell_temperature+=
	      (   theta_temperature)*new_function_values[q_point]+
	      (1.-theta_temperature)*new_function_values[q_point];
	  }
	average_cell_temperature/=n_q_points;
	
	double cell_thermal_conductivity          = -1.E10;
	double cell_total_volumetric_heat_capacity= -1.E10;
	double cell_ice_saturation                = -1.E10;
	material_data(cell->center()[0],average_cell_temperature,
		      cell_thermal_conductivity,cell_total_volumetric_heat_capacity,
		      cell_ice_saturation);
	
	ice_saturation_int.push_back(cell_ice_saturation);
      }
    const Vector<double> ice_saturation(ice_saturation_int.begin(),ice_saturation_int.end());

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution,"solution");
    data_out.add_data_vector(ice_saturation,"ice_saturation");
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
  }

  template <int dim>
  void Heat_Pipe<dim>::update_met_data ()
  {
    if (parameters.boundary_condition_top.compare("first")==0 ||
	parameters.boundary_condition_top.compare("third")==0)
      {
	/*
	 * Originally this function read a file with date (dd/mm/yyyy) and met data
	 * (air temperature, solar radiation, wind speed, etc) (hence the names).
	 * For the purpose of analysing a simplified composite region in 1D with fixed
	 * top boundary condition I'm assuming that we are reading a file with a single
	 * column corresponding to surface temperature at every time step. We can come
	 * back to more complex met data files later on
	 */
	if (met_data.size()==0)
	  {
	    DataTools data_tools;
	    std::vector<std::string> filenames;
	    filenames.push_back(parameters.top_fixed_value_file);
	    data_tools.read_data (filenames,
				  met_data);

	    std::cout << "\tAvailable surface data lines: " << met_data.size()
		      << std::endl << std::endl;

	    std::vector< std::pair<double,double> > table_temperature_surface;
	    std::vector< std::pair<double,double> > table_temperature_room;
	    for (unsigned int i=0; i<met_data.size(); i++)
	      {
		table_temperature_surface.push_back(std::make_pair(met_data[i][0],met_data[i][1]));
		table_temperature_room.push_back(std::make_pair(met_data[i][0],met_data[i][2]));
	      }

	    for (double t=table_temperature_surface[0].first;
		 t<table_temperature_surface[table_temperature_surface.size()-1].first; t+=time_step)
	      {
		std::vector<double> row_interpolated_temperature_surface;
		row_interpolated_temperature_surface
		  .push_back(t);
		row_interpolated_temperature_surface
		  .push_back(data_tools.interpolate_data(table_temperature_surface, t));

		interpolated_temperature_surface.push_back(row_interpolated_temperature_surface);

		std::vector<double> row_interpolated_temperature_room;
		row_interpolated_temperature_room.push_back(t);
		row_interpolated_temperature_room.push_back(data_tools.interpolate_data(table_temperature_room, t));

		interpolated_temperature_room.push_back(row_interpolated_temperature_room);
	      }
	  }
	
	// old_room_temperature    = interpolated_temperature_room[timestep_number-1][1];
	// new_room_temperature    = interpolated_temperature_room[timestep_number  ][1];
	// old_surface_temperature = interpolated_temperature_surface[timestep_number-1][1];
	// new_surface_temperature = interpolated_temperature_surface[timestep_number  ][1];

	double phase=0;
	double average=5;
	double amplitude=2.;
	double period=24.*3600.;
	old_room_temperature    = average+amplitude*cos((2.*M_PI/period)*((timestep_number-1)*time_step-phase));
	new_room_temperature    = average+amplitude*cos((2.*M_PI/period)*( timestep_number   *time_step-phase));
	old_surface_temperature = average+amplitude*cos((2.*M_PI/period)*((timestep_number-1)*time_step-phase));
	new_surface_temperature = average+amplitude*cos((2.*M_PI/period)*( timestep_number   *time_step-phase));
      }
    
    if (parameters.point_source==true)
      {
    	if (point_source_magnitudes.size()==0)
    	  {
    	    DataTools data_tools;
    	    std::vector<std::string> filenames;
    	    filenames.push_back(parameters.point_source_file);
    	    data_tools.read_data (filenames,
    				  point_source_magnitudes);

    	    std::cout << "\n\tPoint source active at: " << parameters.point_source_depth
    		      << "\n\tAvailable point source entries: "
    		      << point_source_magnitudes.size()
    		      << std::endl << std::endl;
    	  }
    	// old_point_source_magnitude =point_source_magnitudes[timestep_number-1][1];
    	// new_point_source_magnitude =point_source_magnitudes[timestep_number  ][1];
    	old_point_source_magnitude=
    	  -point_source_magnitudes[timestep_number-1][1]*sin((2.*M_PI/86400)*((timestep_number-1)*time_step-54000));
    	new_point_source_magnitude=
    	  -point_source_magnitudes[timestep_number  ][1]*sin((2.*M_PI/86400)*((timestep_number  )*time_step-54000));
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
    data_tools.read_data(filenames,
			 initial_condition);
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
      double k=0.;
      double Cp=0.;
      double Si=0.;
      material_data(-parameters.material_0_depth-0.5*parameters.material_0_thickness,
		    25.,k,Cp,Si);
	    
      std::cout.setf( std::ios::fixed);
      std::cout.precision(3);
      std::cout << "\tPosition of material layers:\n"
		<< "\t\tLayer 1: from "
		<< parameters.material_0_depth << " to "
		<< parameters.material_0_depth+parameters.material_0_thickness << "\t"
		<< "k(@25C) :" << k << " W/mK\t"
		<< "Cp(@25C):" << Cp/1.E6 << " MJ/m3K\n";

      material_data(-parameters.material_1_depth-0.5*parameters.material_1_thickness,
		    25.,k,Cp,Si);
      std::cout << "\t\tLayer 2: from "
		<< parameters.material_1_depth << " to "
		<< parameters.material_1_depth+parameters.material_1_thickness << "\t"
		<< "k(@25C) :" << k << " W/mK\t"
		<< "Cp(@25C):" << Cp/1.E6 << " MJ/m3K\n";
	    
      material_data(-parameters.material_2_depth-0.5*parameters.material_2_thickness,
		    25.,k,Cp,Si);
      std::cout << "\t\tLayer 3: from "
		<< parameters.material_2_depth << " to "
		<< parameters.material_2_depth+parameters.material_2_thickness << "\t"
		<< "k(@25C) :" << k << " W/mK\t"
		<< "Cp(@25C):" << Cp/1.E6 << "MJ/m3K\n";

      material_data(-parameters.material_3_depth-0.5*parameters.material_3_thickness,
		    25.,k,Cp,Si);
      std::cout << "\t\tLayer 4: from "
		<< parameters.material_3_depth << " to "
		<< parameters.material_3_depth+parameters.material_3_thickness << "\t"
		<< "k(@25C) :" << k << " W/mK\t"
		<< "Cp(@25C):" << Cp/1.E6 << "MJ/m3K\n";

      material_data(-parameters.material_4_depth-0.5*parameters.material_4_thickness,
		    25.,k,Cp,Si);
      std::cout << "\t\tLayer 5: from "
		<< parameters.material_4_depth << " to "
		<< parameters.material_4_depth+parameters.material_4_thickness << "\t"
		<< "k(@25C) :" << k << " W/mK\t"
		<< "Cp(@25C):" << Cp/1.E6 << "MJ/m3K\n";
    }
    int output_count=0;
    for (timestep_number=1;//,time=time_step;
	 timestep_number<=timestep_number_max;//time<=time_max;
	 ++timestep_number)//time+=time_step
      {
	update_met_data();
	
	int iteration=0;
	double total_error =1.E10;
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
	  }while (std::fabs(total_error)>5E-4);
	
	time+=time_step;
	
	if (parameters.output_data_in_terminal==true)
	  std::cout << "Time step " << timestep_number << "\ttime: " << time/60 << " min\tDt: "
		    << time_step << " s\t#it: " << iteration
		    << "\n";
	
	if (parameters.output_frequency!=0 &&
	    time>output_count*parameters.output_frequency)
	  {
	    output_results();
	    output_count++;
	  }
	fill_output_vectors();
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
