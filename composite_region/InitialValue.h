/***************************************************************************************************
                            Initial Temperature Profile
 ***************************************************************************************************/

template <int dim>
class InitialValue : public Function<dim>
{
public:
  InitialValue(const double constant_temperature_); // constant initial temperature
  // InitialValue(const std::string surface_type_,     // analytical inital temperature
  // 	       const double time_);
  InitialValue(std::vector<std::pair<double,double> > initial_sensor_temperatures_); // experimental initial temperature
  
  virtual double value(const Point<dim> &p,
		       const unsigned int component=0) const;
  virtual void vector_value(const Point<dim> &p,
			    Vector<double> &v) const;
private:
  unsigned int initial_case;
  std::string surface_type;
  double time;
  double constant_temperature;
  std::vector<std::pair<double,double> > initial_sensor_temperatures;
};

template <int dim>
InitialValue<dim>::InitialValue(const double constant_temperature_){
  initial_case=0;
  surface_type="";
  time=0;
  constant_temperature=constant_temperature_;
}

// template <int dim>
// InitialValue<dim>::InitialValue(const std::string surface_type_,
// 			   const double time_){
//   initial_case=1;
//   surface_type=surface_type_;
//   time=time_;
//   constant_temperature=0;
// }

template <int dim>
InitialValue<dim>::InitialValue(std::vector<std::pair<double,double> > initial_sensor_temperatures_){
  initial_case=2;
  surface_type="";
  time=0;
  constant_temperature=0;
  initial_sensor_temperatures=initial_sensor_temperatures_;
}

template <int dim>
double InitialValue<dim>::value (const Point<dim>  &p,
				 const unsigned int /*component*/) const
{
  double val=0;
  
  switch (initial_case)
    {
    case 0:
      val=constant_temperature;
      break;
    case 1:
      {
	// const double BaseTime=400.*365.25*24.*3600.;
	// AnalyticSolution analytic_solution(1.2,840,1960,surface_type);
	// val=analytic_solution.get_value (-1.*p[dim-1],BaseTime + time,5000);
	break;
      }
    case 2:
      {
	double x = 0.0;
	double z = 0.0;
  
	if (dim!=1)
	  {
	    x=p[0];
	    z=p[dim-1];
	  
	    if ((x<=-4.0) && (x>-7.8196))
	      z = z - (x+4.0)*0.7;
	    else if (x<-7.8196)
	      z = z - (-7.8196+4.0)*0.7;
	    else
	      z = z;
	  
	    if (-1.*z > initial_sensor_temperatures.back().first)
	      z=-1.*initial_sensor_temperatures.back().first;
	    val=interpolate_data(initial_sensor_temperatures,-1.*z);
	  }
	else
	  {
	    val=interpolate_data(initial_sensor_temperatures,-1.*p[dim-1]);
	    //std::cout << -1.*p[dim-1] << "\t" << val << std::endl;
	  }
	break;
      }
    default:
      val=-1.e6;
    }
  
  return (val);
}

template <int dim>
void InitialValue<dim>::vector_value (const Point<dim> &p,
				      Vector<double> &v) const
{
  const unsigned int n = this->n_components;
  AssertDimension (v.size(), n);
  
  for (unsigned int i=0;i<this->n_components;++i)
    v(i) = value(p, i); 
}
