namespace Parameters
{  
  using namespace dealii;

  template <int dim>
    struct AllParameters
    {
      AllParameters ();

      unsigned int timestep_number_max;
      double time_step;
      double theta;
      double domain_size;
      double point_source_depth;
      unsigned int number_of_layers;
      unsigned int refinement_level;
      unsigned int output_frequency;

      double material_0_thermal_conductivity_solids;
      double material_1_thermal_conductivity_solids;
      double material_2_thermal_conductivity_solids;
      double material_3_thermal_conductivity_solids;
      double material_4_thermal_conductivity_solids;

      double material_0_density;
      double material_1_density;
      double material_2_density;
      double material_3_density;
      double material_4_density;

      double material_0_specific_heat_capacity;
      double material_1_specific_heat_capacity;
      double material_2_specific_heat_capacity;
      double material_3_specific_heat_capacity;
      double material_4_specific_heat_capacity;

      double material_0_thickness;
      double material_1_thickness;
      double material_2_thickness;
      double material_3_thickness;
      double material_4_thickness;

      double material_0_depth;
      double material_1_depth;
      double material_2_depth;
      double material_3_depth;
      double material_4_depth;

      double material_0_degree_of_saturation;
      double material_1_degree_of_saturation;
      double material_2_degree_of_saturation;
      double material_3_degree_of_saturation;
      double material_4_degree_of_saturation;

      double material_0_porosity;
      double material_1_porosity;
      double material_2_porosity;
      double material_3_porosity;
      double material_4_porosity;

      std::string material_0_thermal_conductivity_relationship;
      std::string material_1_thermal_conductivity_relationship;
      std::string material_2_thermal_conductivity_relationship;
      std::string material_3_thermal_conductivity_relationship;
      std::string material_4_thermal_conductivity_relationship;

      std::string material_0_name;
      std::string material_1_name;
      std::string material_2_name;
      std::string material_3_name;
      std::string material_4_name;
      
      double freezing_point;
      double alpha;
      double latent_heat;
      double reference_temperature;
      double heat_loss_factor;

      double density_ice;
      double density_air;
      double density_liquids;

      double specific_heat_capacity_ice;
      double specific_heat_capacity_air;
      double specific_heat_capacity_liquids;

      double thermal_conductivity_liquids;
      double thermal_conductivity_air;

      bool fixed_at_bottom;
      double bottom_fixed_value;
      bool fixed_at_top;
      bool point_source;
      bool output_data_in_terminal;

      std::string boundary_condition_top;

      std::string top_fixed_value_file;
      std::string initial_condition_file;
      std::string depths_file;
      std::string point_source_file;
      std::string output_directory;
      std::string output_file;

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };

  template <int dim>
    AllParameters<dim>::AllParameters ()
    {
      timestep_number_max=0;
      time_step=0.;
      theta=0.;
      domain_size=0.;
      point_source_depth=0.;
      number_of_layers=0;
      refinement_level=0;
      output_frequency=0;

      material_0_thermal_conductivity_solids=0.;
      material_1_thermal_conductivity_solids=0.;
      material_2_thermal_conductivity_solids=0.;
      material_3_thermal_conductivity_solids=0.;
      material_4_thermal_conductivity_solids=0.;

      material_0_density=0.;
      material_1_density=0.;
      material_2_density=0.;
      material_3_density=0.;
      material_4_density=0.;

      material_0_specific_heat_capacity=0.;
      material_1_specific_heat_capacity=0.;
      material_2_specific_heat_capacity=0.;
      material_3_specific_heat_capacity=0.;
      material_4_specific_heat_capacity=0.;

      material_0_thickness=0.;
      material_1_thickness=0.;
      material_2_thickness=0.;
      material_3_thickness=0.;
      material_4_thickness=0.;

      material_0_depth=0.;
      material_1_depth=0.;
      material_2_depth=0.;
      material_3_depth=0.;
      material_4_depth=0.;

      material_0_degree_of_saturation=0.;
      material_1_degree_of_saturation=0.;
      material_2_degree_of_saturation=0.;
      material_3_degree_of_saturation=0.;
      material_4_degree_of_saturation=0.;

      material_0_porosity=0.;
      material_1_porosity=0.;
      material_2_porosity=0.;
      material_3_porosity=0.;
      material_4_porosity=0.;

      freezing_point=0.;
      alpha=0.;
      latent_heat=0.;
      reference_temperature=0.;
      heat_loss_factor=0.;

      density_ice=0.;
      density_air=0.;
      density_liquids=0.;

      specific_heat_capacity_ice=0.;
      specific_heat_capacity_air=0.;
      specific_heat_capacity_liquids=0.;

      thermal_conductivity_liquids=0.;
      thermal_conductivity_air=0.;

      fixed_at_bottom=false;
      bottom_fixed_value=0.;
      fixed_at_top=false;
      point_source=false;
      output_data_in_terminal=false;
    }

  template <int dim>
    void
    AllParameters<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("time stepping");
      {
	prm.declare_entry("time step", "3600",
			  Patterns::Double(0),
			  "simulation time step");
	prm.declare_entry("timestep number max", "70055",
			  Patterns::Integer(),
			  "max number of timesteps to execute");
	prm.declare_entry("theta scheme value", "0.5",
			  Patterns::Double(0,1),
			  "value for theta that interpolated between explicit "
			  "Euler (theta=0), Crank-Nicolson (theta=0.5), and "
			  "implicit Euler (theta=1).");
      }
      prm.leave_subsection();

      prm.enter_subsection("geometric data");
      {
	prm.declare_entry("domain size", "20",
			  Patterns::Double(0),
			  "size of domain in m");
	prm.declare_entry("number of layers", "1",
			  Patterns::Integer(1,5),
			  "number of layers composing the"
			  " domain");
	prm.declare_entry("material 0 thickness", "0.20",
			  Patterns::Double(0),
			  "thickness of material 0 in m");
	prm.declare_entry("material 0 depth", "1.00",
			  Patterns::Double(0),
			  "depth of material 0 in m");
	prm.declare_entry("material 1 thickness", "0.20",
			  Patterns::Double(0),
			  "thickness of material 1 in m");
	prm.declare_entry("material 1 depth", "1.00",
			  Patterns::Double(0),
			  "depth of material 1 in m");
	prm.declare_entry("material 2 thickness", "0.20",
			  Patterns::Double(0),
			  "thickness of material 2 in m");
	prm.declare_entry("material 2 depth", "1.00",
			  Patterns::Double(0),
			  "depth of material 2 in m");
	prm.declare_entry("material 3 thickness", "0.20",
			  Patterns::Double(0),
			  "thickness of material 3 in m");
	prm.declare_entry("material 3 depth", "1.00",
			  Patterns::Double(0),
			  "depth of material 3 in m");
	prm.declare_entry("material 4 thickness", "0.20",
			  Patterns::Double(0),
			  "thickness of material 4 in m");
	prm.declare_entry("material 4 depth", "1.00",
			  Patterns::Double(0),
			  "depth of material 4 in m");
	prm.declare_entry("refinement level", "5",
			  Patterns::Integer(),
			  "number of cells as in 2^n");
      }
      prm.leave_subsection();


      prm.enter_subsection("material data");
      {
	prm.declare_entry("air density",
			  "0.",Patterns::Double(0),
			  "density of air in kg/m3");
	prm.declare_entry("liquids density",
			  "0.",Patterns::Double(0),
			  "density of liquids in kg/m3");
	prm.declare_entry("air specific heat capacity",
			  "0.",Patterns::Double(0),
			  "specific capacity of air in J/kgK");
	prm.declare_entry("liquids specific heat capacity",
			  "0.",Patterns::Double(0),
			  "specific capacity of liquids in J/kgK");
	prm.declare_entry("ice density",
			  "0.",Patterns::Double(0),
			  "density of ice in kg/m3");
	prm.declare_entry("ice specific heat capacity",
			  "0.",Patterns::Double(0),
			  "specific capacity of ice in J/kgK");
	prm.declare_entry("initial temperature",
			  "0.",Patterns::Double(0),
			  "laboratory temperature");
	prm.declare_entry("reference temperature",
			  "0.",Patterns::Double(),
			  "reference temperature");
	prm.declare_entry("freezing point",
			  "0.",Patterns::Double(-20.),
			  "freezing point of pore water");
	prm.declare_entry("alpha",
			  "0.",Patterns::Double(-10.,0),
			  "alpha of pore water");
	prm.declare_entry("latent heat",
			  "0.",Patterns::Double(0),
			  "latent heat of fusion");
	prm.declare_entry("liquids thermal conductivity",
			  "0.",Patterns::Double(0),
			  "thermal conductivity of liquids in W/mK");
	prm.declare_entry("air thermal conductivity",
			  "0.",Patterns::Double(0),
			  "thermal conductivity of air in W/mK");
	/*
	 * layer 0
	 * */
	prm.declare_entry("material 0 degree of saturation",
			  "0.",Patterns::Double(0.,1.),
			  "degree of saturation of soil layer 0");
	prm.declare_entry("material 0 porosity",
			  "0.",Patterns::Double(0.,1.),
			  "porosity of soil layer 0");
	prm.declare_entry("material 0 thermal conductivity solids",
			  "0.",Patterns::Double(0.),
			  "thermal conductivity of material 0 in W/mK");
	prm.declare_entry("material 0 density",
			  "0.",Patterns::Double(0.),
			  "density of soil layer 0 in kg/m3");
	prm.declare_entry("material 0 specific heat capacity",
			  "0.",Patterns::Double(0.),
			  "specific capacity of soil layer 0 in J/kgK");
	prm.declare_entry("material 0 thermal conductivity relationship",
			  "",Patterns::Anything(),
			  "string defining the theoretical relationship "
			  "to estimate the thermal conductivity of the layer. "
			  "Three expressions are currently defined based on "
			  "the work of 'hugh' (2012); 'donazzi' (1979); and "
			  "a 'bulk' relation that uses the value provided to"
			  "thermal_conductivity_solids as it is.");
	prm.declare_entry("material 0 name",
			  "",Patterns::Anything(),
			  "Name of the material comprising material 0. The "
			  "name is searched in a Map and if found the "
			  "corresponding properties of the solid particles "
			  "are accessed.");
	/*
	 * layer 1
	 * */
	prm.declare_entry("material 1 degree of saturation",
			  "0.",Patterns::Double(0.,1.),
			  "degree of saturation of soil layer 1");
	prm.declare_entry("material 1 porosity",
			  "0.",Patterns::Double(0.,1.),
			  "porosity of soil layer 1");
	prm.declare_entry("material 1 thermal conductivity solids",
			  "0.",Patterns::Double(0.),
			  "thermal conductivity of soil layer 1 in W/mK");
	prm.declare_entry("material 1 density",
			  "0.",Patterns::Double(0.),
			  "density of soil layer 1 in kg/m3");
	prm.declare_entry("material 1 specific heat capacity",
			  "0.",Patterns::Double(0.),
			  "specific capacity of soil layer 1 in J/kgK");
	prm.declare_entry("material 1 thermal conductivity relationship",
			  "",Patterns::Anything(),
			  "string defining the theoretical relationship "
			  "to estimate the thermal conductivity of the layer. "
			  "Three expressions are currently defined based on "
			  "the work of 'hugh' (2012); 'donazzi' (1979); and "
			  "a 'bulk' relation that uses the value provided to"
			  "thermal_conductivity_solids as it is.");
	prm.declare_entry("material 1 name",
			  "",Patterns::Anything(),
			  "Name of the material comprising material 1. The "
			  "name is searched in a Map and if found the "
			  "corresponding properties of the solid particles "
			  "are accessed.");
	/*
	 * layer 2
	 * */
	prm.declare_entry("material 2 degree of saturation",
			  "0.",Patterns::Double(0.,1.),
			  "degree of saturation of soil layer 2");
	prm.declare_entry("material 2 porosity",
			  "0.",Patterns::Double(0.,1.),
			  "porosity of soil layer 2");
	prm.declare_entry("material 2 thermal conductivity solids",
			  "0.",Patterns::Double(0.),
			  "thermal conductivity of soil layer 2 in W/mK");
	prm.declare_entry("material 2 density",
			  "0.",Patterns::Double(0.),
			  "density of soil layer 2 in kg/m3");
	prm.declare_entry("material 2 specific heat capacity",
			  "0.",Patterns::Double(0.),
			  "specific capacity of soil layer 2 in J/kgK");
	prm.declare_entry("material 2 thermal conductivity relationship",
			  "",Patterns::Anything(),
			  "string defining the theoretical relationship "
			  "to estimate the thermal conductivity of the layer. "
			  "Three expressions are currently defined based on "
			  "the work of 'hugh' (2012); 'donazzi' (1979); and "
			  "a 'bulk' relation that uses the value provided to"
			  "thermal_conductivity_solids as it is.");
	prm.declare_entry("material 2 name",
			  "",Patterns::Anything(),
			  "Name of the material comprising material 2. The "
			  "name is searched in a Map and if found the "
			  "corresponding properties of the solid particles "
			  "are accessed.");
	/*
	 * layer 3
	 * */
	prm.declare_entry("material 3 degree of saturation",
			  "0.",Patterns::Double(0.,1.),
			  "degree of saturation of soil layer 3");
	prm.declare_entry("material 3 porosity",
			  "0.",Patterns::Double(0.,1.),
			  "porosity of soil layer 3");
	prm.declare_entry("material 3 thermal conductivity solids",
			  "0.",Patterns::Double(0.),
			  "thermal conductivity of soil layer 3 in W/mK");
	prm.declare_entry("material 3 density",
			  "0.",Patterns::Double(0.),
			  "density of soil layer 3 in kg/m3");
	prm.declare_entry("material 3 specific heat capacity",
			  "0.",Patterns::Double(0.),
			  "specific capacity of soil layer 3 in J/kgK");
	prm.declare_entry("material 3 thermal conductivity relationship",
			  "",Patterns::Anything(),
			  "string defining the theoretical relationship "
			  "to estimate the thermal conductivity of the layer. "
			  "Three expressions are currently defined based on "
			  "the work of 'hugh' (2012); 'donazzi' (1979); and "
			  "a 'bulk' relation that uses the value provided to"
			  "thermal_conductivity_solids as it is.");
	prm.declare_entry("material 3 name",
			  "",Patterns::Anything(),
			  "Name of the material comprising material 3. The "
			  "name is searched in a Map and if found the "
			  "corresponding properties of the solid particles "
			  "are accessed.");
	/*
	 * layer 4
	 * */
	prm.declare_entry("material 4 degree of saturation",
			  "0.",Patterns::Double(0.,1.),
			  "degree of saturation of soil layer 4");
	prm.declare_entry("material 4 porosity",
			  "0.",Patterns::Double(0.,1.),
			  "porosity of soil layer 4");
	prm.declare_entry("material 4 thermal conductivity solids",
			  "0.",Patterns::Double(0.),
			  "thermal conductivity of soil layer 4 in W/mK");
	prm.declare_entry("material 4 density",
			  "0.",Patterns::Double(0.),
			  "density of soil layer 4 in kg/m3");
	prm.declare_entry("material 4 specific heat capacity",
			  "0.",Patterns::Double(0.),
			  "specific capacity of soil layer 4 in J/kgK");
	prm.declare_entry("material 4 thermal conductivity relationship",
			  "",Patterns::Anything(),
			  "string defining the theoretical relationship "
			  "to estimate the thermal conductivity of the layer. "
			  "Three expressions are currently defined based on "
			  "the work of 'hugh' (2012); 'donazzi' (1979); and "
			  "a 'bulk' relation that uses the value provided to"
			  "thermal_conductivity_solids as it is.");
	prm.declare_entry("material 4 name",
			  "",Patterns::Anything(),
			  "Name of the material comprising material 4. The "
			  "name is searched in a Map and if found the "
			  "corresponding properties of the solid particles "
			  "are accessed.");
      }
      prm.leave_subsection();


      prm.enter_subsection("boundary conditions");
      {
	prm.declare_entry("fixed at bottom", "false",
			  Patterns::Bool(),
			  "if true, the bottom boundary condtion is fixed");
	prm.declare_entry("bottom fixed value", "10.",
			  Patterns::Double(),
			  "value at bottom boundary for fixed conditions");
	prm.declare_entry("fixed at top", "true",
			  Patterns::Bool(),
			  "if true, the top boundary condtion is fixed");
	prm.declare_entry("point source", "false",
			  Patterns::Bool(),
			  "if true, a point source is present in the domain");
	prm.declare_entry("top fixed value file","surface_temperature.txt",
			  Patterns::Anything(),
			  "file containing values for the top boundary"
			  "in case of fixed conditions.");
	prm.declare_entry("initial condition file","initial_condition.txt",
			  Patterns::Anything(),
			  "file containing values of temperature to be "
			  "used as initial condition.");
	prm.declare_entry("depths file","depths.txt",
			  Patterns::Anything(),
			  "file containing coordinates where data will "
			  "be extracted.");
	prm.declare_entry("point source file","point_source_file.txt",
			  Patterns::Anything(),
			  "file containing values for the magnitude "
			  "of a point source in the domain.");
	prm.declare_entry("point source depth", "7.00",
			  Patterns::Double(0),
			  "depth of point souce in m");
	prm.declare_entry("heat loss factor","0.0",
			  Patterns::Double(0.),
			  "Heat loss factor in W/m3K. This is estimated integrating a"
			  " surface heat transfer coefficient (W/m2K) over the boundary"
			  " and dividing it by the volume of the domain");
	prm.declare_entry("boundary condition top","first",
			  Patterns::Anything(),
			  "set to first, second or third to set the "
			  "corresponding boundary condition at the top");
      }
      prm.leave_subsection();

      prm.enter_subsection("other options");
      {
	prm.declare_entry("output frequency", "1",
			  Patterns::Integer(),"number that defines the"
			  "number of time steps between outputs");
	prm.declare_entry("output directory", "output",
			  Patterns::Anything(),
			  "name of output directory to store all"
			  "output files.");
	prm.declare_entry("output file", "output_data.txt",
			  Patterns::Anything(), "Defines the name of the filename "
			  "to store the temperatures at the points defined "
			  "in the file 'depths_file'");
	prm.declare_entry("output data in terminal", "true",
			  Patterns::Bool(),"if true, the program will generate output "
			  "in the terminal. Set to false to avoid cluttering "
			  "and speed up a bit the program.");
      }
      prm.leave_subsection();
    }

  template <int dim>
    void AllParameters<dim>::parse_parameters (ParameterHandler &prm)
    {
      //mesh_filename = prm.get("mesh");

      prm.enter_subsection("time stepping");
      {
	time_step           = prm.get_double ("time step");
	timestep_number_max = prm.get_integer("timestep number max");
	theta               = prm.get_double ("theta scheme value");
      }
      prm.leave_subsection();


      prm.enter_subsection("geometric data");
      {
	domain_size           = prm.get_double ("domain size");
	number_of_layers      = prm.get_integer("number of layers");
	refinement_level      = prm.get_integer("refinement level");
	material_0_depth      = prm.get_double ("material 0 depth");
	material_0_thickness  = prm.get_double ("material 0 thickness");
	material_1_depth      = prm.get_double ("material 1 depth");
	material_1_thickness  = prm.get_double ("material 1 thickness");
	material_2_depth      = prm.get_double ("material 2 depth");
	material_2_thickness  = prm.get_double ("material 2 thickness");
	material_3_depth      = prm.get_double ("material 3 depth");
	material_3_thickness  = prm.get_double ("material 3 thickness");
	material_4_depth      = prm.get_double ("material 4 depth");
	material_4_thickness  = prm.get_double ("material 4 thickness");
      }
      prm.leave_subsection();


      prm.enter_subsection("material data");
      {
	density_air                       = prm.get_double ("air density");
	density_liquids                   = prm.get_double ("liquids density");
	specific_heat_capacity_air        = prm.get_double ("air specific heat capacity");
	specific_heat_capacity_liquids    = prm.get_double ("liquids specific heat capacity");
	specific_heat_capacity_ice        = prm.get_double ("ice specific heat capacity");
	density_ice                       = prm.get_double ("ice density");
	reference_temperature             = prm.get_double ("reference temperature");
	freezing_point                    = prm.get_double ("freezing point");
	alpha                             = prm.get_double ("alpha");
	latent_heat                       = prm.get_double ("latent heat");
	thermal_conductivity_liquids      = prm.get_double ("liquids thermal conductivity");
	thermal_conductivity_air          = prm.get_double ("air thermal conductivity");
	material_0_degree_of_saturation   = prm.get_double ("material 0 degree of saturation");
	material_0_porosity               = prm.get_double ("material 0 porosity");
	material_0_density                = prm.get_double ("material 0 density");
	material_0_specific_heat_capacity = prm.get_double ("material 0 specific heat capacity");
	material_0_thermal_conductivity_solids
	  =prm.get_double ("material 0 thermal conductivity solids");
	material_0_thermal_conductivity_relationship
	  =prm.get("material 0 thermal conductivity relationship");
	material_0_name                   = prm.get("material 0 name");
	  
	material_1_degree_of_saturation   = prm.get_double ("material 1 degree of saturation");
	material_1_porosity               = prm.get_double ("material 1 porosity");
	material_1_density                = prm.get_double ("material 1 density");
	material_1_specific_heat_capacity = prm.get_double ("material 1 specific heat capacity");
	material_1_thermal_conductivity_solids
	  =prm.get_double ("material 1 thermal conductivity solids");
	material_1_thermal_conductivity_relationship
	  =prm.get("material 1 thermal conductivity relationship");
	material_1_name                   = prm.get("material 1 name");
	
	material_2_degree_of_saturation   = prm.get_double ("material 2 degree of saturation");
	material_2_porosity               = prm.get_double ("material 2 porosity");
	material_2_density                = prm.get_double ("material 2 density");
	material_2_specific_heat_capacity = prm.get_double ("material 2 specific heat capacity");
	material_2_thermal_conductivity_solids
	  = prm.get_double ("material 2 thermal conductivity solids");
	material_2_thermal_conductivity_relationship
	  =prm.get("material 2 thermal conductivity relationship");
	material_2_name                   = prm.get("material 2 name");
	
	material_3_degree_of_saturation   = prm.get_double ("material 3 degree of saturation");
	material_3_porosity               = prm.get_double ("material 3 porosity");
	material_3_density                = prm.get_double ("material 3 density");
	material_3_specific_heat_capacity = prm.get_double ("material 3 specific heat capacity");
	material_3_thermal_conductivity_solids
	  =prm.get_double ("material 3 thermal conductivity solids");
	material_3_thermal_conductivity_relationship
	  =prm.get("material 3 thermal conductivity relationship");
	material_3_name                   = prm.get("material 3 name");
	
	material_4_degree_of_saturation   = prm.get_double ("material 4 degree of saturation");
	material_4_porosity               = prm.get_double ("material 4 porosity");
	material_4_density                = prm.get_double ("material 4 density");
	material_4_specific_heat_capacity = prm.get_double ("material 4 specific heat capacity");
	material_4_thermal_conductivity_solids
	  =prm.get_double ("material 4 thermal conductivity solids");
	material_4_thermal_conductivity_relationship
	  =prm.get("material 4 thermal conductivity relationship");
	material_4_name                   = prm.get("material 4 name");
      }
      prm.leave_subsection();


      prm.enter_subsection("boundary conditions");
      {
	fixed_at_bottom           = prm.get_bool   ("fixed at bottom");
	bottom_fixed_value        = prm.get_double ("bottom fixed value");
	fixed_at_top              = prm.get_bool   ("fixed at top");
	point_source              = prm.get_bool   ("point source");
	point_source_depth        = prm.get_double ("point source depth");
	top_fixed_value_file      = prm.get        ("top fixed value file"); 
	initial_condition_file    = prm.get        ("initial condition file");
	depths_file               = prm.get        ("depths file");
	point_source_file         = prm.get        ("point source file");
	heat_loss_factor          = prm.get_double ("heat loss factor");
	boundary_condition_top    = prm.get        ("boundary condition top");
      }
      prm.leave_subsection();

      prm.enter_subsection("other options");
      {
	output_frequency	= prm.get_integer("output frequency");
	output_directory	= prm.get	 ("output directory");
	output_file         = prm.get    ("output file");
	output_data_in_terminal=prm.get_bool("output data in terminal");
      }
      prm.leave_subsection();
    }
}
