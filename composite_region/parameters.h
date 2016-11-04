
namespace Parameters
{  
template <int dim>
struct AllParameters
{
	AllParameters ();

	unsigned int timestep_number_max;
	double time_step;
	double theta;
	double domain_size;
	double point_source_depth;
	unsigned int refinement_level;
	unsigned int output_frequency;

	double material_0_thermal_conductivity;
	double material_1_thermal_conductivity;
	double material_2_thermal_conductivity;
	double material_3_thermal_conductivity;
	double material_4_thermal_conductivity;

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

	double freezing_point;
	double alpha;
	double latent_heat;
	double initial_temperature;
	double reference_temperature;

	double density_ice;
	double density_air;
	double density_liquids;

	double specific_heat_capacity_ice;
	double specific_heat_capacity_air;
	double specific_heat_capacity_liquids;

	bool fixed_at_bottom;
	double bottom_fixed_value;
	bool fixed_at_top;
	bool point_source;

	std::string top_fixed_value_file;
	std::string initial_condition_file;
	std::string depths_file;
	std::string point_source_file;
	std::string output_directory;

	static void declare_parameters (ParameterHandler &prm);
	void parse_parameters (ParameterHandler &prm);
};

template <int dim>
AllParameters<dim>::AllParameters ()
{}

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
		prm.declare_entry("point source depth", "7.00",
				Patterns::Double(0),
				"depth of point souce in m");
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
				"specific capacity of air in J/mK");
		prm.declare_entry("liquids specific heat capacity",
				"0.",Patterns::Double(0),
				"specific capacity of liquids in J/mK");
		prm.declare_entry("ice density",
				"0.",Patterns::Double(0),
				"density of ice in kg/m3");
		prm.declare_entry("ice specific heat capacity",
				"0.",Patterns::Double(0),
				"specific capacity of ice in J/mK");
		prm.declare_entry("initial temperature",
				"0.",Patterns::Double(0),
				"laboratory temperature");
		prm.declare_entry("reference temperature",
				"0.",Patterns::Double(0),
				"reference temperature");
		prm.declare_entry("freezing point",
				"0.",Patterns::Double(0),
				"freezing point of pore water");
		prm.declare_entry("alpha",
				"0.",Patterns::Double(-10.,0),
				"alpha of pore water");
		prm.declare_entry("latent heat",
				"0.",Patterns::Double(0),
				"latent heat of fusion");
		/*
		 * layer 0
		 * */
		prm.declare_entry("material 0 degree of saturation",
				"0.",Patterns::Double(0.,1.),
				"degree of saturation of soil layer 0");
		prm.declare_entry("material 0 porosity",
				"0.",Patterns::Double(0.,1.),
				"porosity of soil layer 0");
		prm.declare_entry("material 0 thermal conductivity",
				"0.",Patterns::Double(0.),
				"thermal conductivity of material 0 in W/mK");
		prm.declare_entry("material 0 density",
				"0.",Patterns::Double(0.),
				"density of soil layer 0 in kg/m3");
		prm.declare_entry("material 0 specific heat capacity",
				"0.",Patterns::Double(0.),
				"specific capacity of soil layer 0 in J/kgK");
		/*
		 * layer 1
		 * */
		prm.declare_entry("material 1 degree of saturation",
				"0.",Patterns::Double(0.,1.),
				"degree of saturation of soil layer 1");
		prm.declare_entry("material 1 porosity",
				"0.",Patterns::Double(0.,1.),
				"porosity of soil layer 1");
		prm.declare_entry("material 1 thermal conductivity",
				"0.",Patterns::Double(0.),
				"thermal conductivity of soil layer 1 in W/mK");
		prm.declare_entry("material 1 density",
				"0.",Patterns::Double(0.),
				"density of soil layer 1 in kg/m3");
		prm.declare_entry("material 1 specific heat capacity",
				"0.",Patterns::Double(0.),
				"specific capacity of soil layer 1 in J/kgK");
		/*
		 * layer 2
		 * */
		prm.declare_entry("material 2 degree of saturation",
				"0.",Patterns::Double(0.,1.),
				"degree of saturation of soil layer 2");
		prm.declare_entry("material 2 porosity",
				"0.",Patterns::Double(0.,1.),
				"porosity of soil layer 2");
		prm.declare_entry("material 2 thermal conductivity",
				"0.",Patterns::Double(0.),
				"thermal conductivity of soil layer 2 in W/mK");
		prm.declare_entry("material 2 density",
				"0.",Patterns::Double(0.),
				"density of soil layer 2 in kg/m3");
		prm.declare_entry("material 2 specific heat capacity",
				"0.",Patterns::Double(0.),
				"specific capacity of soil layer 2 in J/kgK");
		/*
		 * layer 3
		 * */
		prm.declare_entry("material 3 degree of saturation",
				"0.",Patterns::Double(0.,1.),
				"degree of saturation of soil layer 3");
		prm.declare_entry("material 3 porosity",
				"0.",Patterns::Double(0.,1.),
				"porosity of soil layer 3");
		prm.declare_entry("material 3 thermal conductivity",
				"0.",Patterns::Double(0.),
				"thermal conductivity of soil layer 3 in W/mK");
		prm.declare_entry("material 3 density",
				"0.",Patterns::Double(0.),
				"density of soil layer 3 in kg/m3");
		prm.declare_entry("material 3 specific heat capacity",
				"0.",Patterns::Double(0.),
				"specific capacity of soil layer 3 in J/kgK");
		/*
		 * layer 4
		 * */
		prm.declare_entry("material 4 degree of saturation",
				"0.",Patterns::Double(0.,1.),
				"degree of saturation of soil layer 4");
		prm.declare_entry("material 4 porosity",
				"0.",Patterns::Double(0.,1.),
				"porosity of soil layer 4");
		prm.declare_entry("material 4 thermal conductivity",
				"0.",Patterns::Double(0.),
				"thermal conductivity of soil layer 4 in W/mK");
		prm.declare_entry("material 4 density",
				"0.",Patterns::Double(0.),
				"density of soil layer 4 in kg/m3");
		prm.declare_entry("material 4 specific heat capacity",
				"0.",Patterns::Double(0.),
				"specific capacity of soil layer 4 in J/kgK");
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
		point_source_depth    = prm.get_double ("point source depth");
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
		initial_temperature               = prm.get_double ("initial temperature");
		reference_temperature             = prm.get_double ("reference temperature");
		freezing_point                    = prm.get_double ("freezing point");
		alpha                             = prm.get_double ("alpha");
		latent_heat                       = prm.get_double ("latent heat");
		material_0_degree_of_saturation   = prm.get_double ("material 0 degree of saturation");
		material_0_porosity               = prm.get_double ("material 0 porosity");
		material_0_thermal_conductivity   = prm.get_double ("material 0 thermal conductivity");
		material_0_density                = prm.get_double ("material 0 density");
		material_0_specific_heat_capacity = prm.get_double ("material 0 specific heat capacity");
		material_1_degree_of_saturation   = prm.get_double ("material 1 degree of saturation");
		material_1_porosity               = prm.get_double ("material 1 porosity");
		material_1_thermal_conductivity   = prm.get_double ("material 1 thermal conductivity");
		material_1_density                = prm.get_double ("material 1 density");
		material_1_specific_heat_capacity = prm.get_double ("material 1 specific heat capacity");
		material_2_degree_of_saturation   = prm.get_double ("material 2 degree of saturation");
		material_2_porosity               = prm.get_double ("material 2 porosity");
		material_2_thermal_conductivity   = prm.get_double ("material 2 thermal conductivity");
		material_2_density                = prm.get_double ("material 2 density");
		material_2_specific_heat_capacity = prm.get_double ("material 2 specific heat capacity");
		material_3_degree_of_saturation   = prm.get_double ("material 3 degree of saturation");
		material_3_porosity               = prm.get_double ("material 3 porosity");
		material_3_thermal_conductivity   = prm.get_double ("material 3 thermal conductivity");
		material_3_density                = prm.get_double ("material 3 density");
		material_3_specific_heat_capacity = prm.get_double ("material 3 specific heat capacity");
		material_4_degree_of_saturation   = prm.get_double ("material 4 degree of saturation");
		material_4_porosity               = prm.get_double ("material 4 porosity");
		material_4_thermal_conductivity   = prm.get_double ("material 4 thermal conductivity");
		material_4_density                = prm.get_double ("material 4 density");
		material_4_specific_heat_capacity = prm.get_double ("material 4 specific heat capacity");
	}
	prm.leave_subsection();


	prm.enter_subsection("boundary conditions");
	{
		fixed_at_bottom        = prm.get_bool   ("fixed at bottom");
		bottom_fixed_value     = prm.get_double ("bottom fixed value");
		fixed_at_top           = prm.get_bool   ("fixed at top");
		point_source           = prm.get_bool   ("point source");
		top_fixed_value_file   = prm.get        ("top fixed value file");
		initial_condition_file = prm.get        ("initial condition file");
		depths_file            = prm.get        ("depths file");
		point_source_file      = prm.get        ("point source file");
	}
	prm.leave_subsection();

	prm.enter_subsection("other options");
	{
		output_frequency	= prm.get_integer("output frequency");
		output_directory	= prm.get		 ("output directory");
	}
	prm.leave_subsection();
}
}
