void print_data (int lines_to_print,
		 const std::vector< std::vector<int> >    &date_and_time,
		 const std::vector< std::vector<double> > &data,
		 std::ostream &outFile) 
{
  // check if the containers have the same size
  if(data.size() != date_and_time.size())
    {
      for (unsigned int i=0; i<date_and_time[0].size(); i++)
	std::cout << date_and_time[0][i] << "\t";
      for (unsigned int i=0; i<data[0].size(); i++)
	std::cout << data[0][i] << "\t";
      std::cout << std::endl;
      
      std::cout << "Error in print_data function. Vector mismatch." << std::endl;
      std::cout << "Vector 1 size : " << date_and_time.size() << std::endl;
      std::cout << "Vector 2 size : " << data.size() << std::endl;      
      throw -1;
    }
  // If we call the function with 0 as number of lines
  // to print, we will print them all. Not really
  // straightforward
  if (lines_to_print == 0)
    lines_to_print = data.size();
  // we are assuming 6 entries in each line of date and time container,
  // 3 for date and 3 for time. Let's check for this
  // if ((date_and_time[0].size() != 6))
  //   throw 4;
  
  if (lines_to_print<0)
    {
      std::cout << "Error. Wrong number of lines to print" << std::endl;
      throw 100;
    }
  
  // print formated output to selected ostream
  for (unsigned int i = 0; i<(unsigned int)lines_to_print; i++)
    {
      outFile << std::setw(2) << std::setfill('0') << date_and_time[i][0] << "/" 
	      << std::setw(2) << std::setfill('0') << date_and_time[i][1] << "/" 
	      << std::setw(2) << std::setfill('0') << date_and_time[i][2] << "\t"
	      << std::setw(2) << std::setfill('0') << date_and_time[i][3] << ":" 
	      << std::setw(2) << std::setfill('0') << date_and_time[i][4] << ":" 
	      << std::setw(2) << std::setfill('0') << date_and_time[i][5];
      
      outFile.setf( std::ios::fixed, std::ios::floatfield );
      for (unsigned int j = 0; j<data[0].size(); j++)
	{
	  if(fabs(data[i][j])<0.001)
	    outFile << "\t" << std::setw(7) << std::setfill(' ') << std::scientific << data[i][j];
	  else
	    outFile << "\t" << std::setw(7) << std::setfill(' ') << std::setprecision(2) << std::fixed << data[i][j];
	}
      // outFile << "\t";
      // for (unsigned int j = 0; j<met_data[0].size(); j++)
      // 	outFile << "\t" /*<< std::setw(8)*/ << std::setfill(' ') << std::setprecision(3) << met_data[i][j];

      outFile << "\n";
    }
}

// //**************************************************************************************************
// void read_data (const std::vector< std::string >   &filenames,
// 		std::vector< std::vector<int> >    &date_and_time,
// 		std::vector< std::vector<double> > &data,
// 		bool day_number_column = true)  // by default we assume that the file has a column for the day number
// {
//   // We clear the containers that will store met data and date data from each
//   // set of files
//   data.clear();
//   date_and_time.clear();
//   // In this loop, we will go through each set of files, reading them,
//   // storing the entries in temporal vectors and storing those vectors
//   // in the above containers

//   for (unsigned int j=0; j<filenames.size(); j++)
//     {
//       // std::cout << "m2" << std::endl;
//       // Open files and check that were correctly opened
//       std::ifstream file (filenames[j].c_str());
//       if (!file.is_open())
// 	throw 2;

//       // Print the file names we are currently working on
//       // std::cout << " \t Reading input file: "
//       // 		<< filenames[j]
//       // 		<< std::endl;
//       // Define some variables that will help us with the
//       // reading process
//       std::string line_file;
//       std::string file_token;
//       std::stringstream file_iss;
//       // We read each line in each file. If EOF is reach
//       // in any of the files, the loop is finished
//       while (std::getline(file,line_file)) 
// 	{
// 	  std::vector<int>    row_date_and_time;
// 	  std::vector<double> row_data;
	    
// 	  file_iss << line_file;
	    
// 	  std::stringstream file_iss1;
// 	  int c = 0;
	  
// 	  // int year = 0;
// 	  // if (j<=8)
// 	  //   year=0;
// 	  // if (j>8 && j<=16)
// 	  //   year=1;

// 	  while (std::getline(file_iss,file_token,'\t'))
// 	    {
// 	      file_iss1 << file_token;
// 	      if (c==0)
// 		while (std::getline(file_iss1,file_token,'/')) row_date_and_time.push_back(atoi(file_token.c_str()));
// 		// {
// 		  // int cc=0;

// 		    // {

// 		    //   if (cc<2)
// 		    // 	row_date_and_time.push_back(atoi(file_token.c_str()));
// 		    //   else
// 		    // 	row_date_and_time.push_back(atoi(file_token.c_str())+year);

// 		    //   cc++;
// 		    // }
// 		// }
// 	      if (c==1)
// 		while (std::getline(file_iss1,file_token,':')) row_date_and_time.push_back(atoi(file_token.c_str())); 
// 	      if (c>1)
// 		while (std::getline(file_iss1,file_token,':'))
// 		  {
// 		    if (day_number_column && c>2)
// 		      row_data.push_back(atof(file_token.c_str()));
// 		    else if (!day_number_column)
// 		      row_data.push_back(atof(file_token.c_str()));
// 		  }
// 	      file_iss1.clear();
// 	      c++;
// 	    }
// 	  file_iss.clear();
// 	  // Store temporal vectors to corresponding containers
// 	  date_and_time.push_back (row_date_and_time);
// 	  data.push_back          (row_data);
// 	}
//       // We finished reading a pair of files. Let's check if we actually reach EOF in both files.
//       // If not, that means that they have a different number of entries. That's not good of course.
//       // if (line_file.bad()) 
//       //   throw 6;
//       // else if (!line_file.eof())
//       //   throw 7;
//       // else 
//       //   {
//       //     // format error (not possible with getline but possible with operator>>)
//       //     // or end of file (can't make the difference)
//       //   }
//       // close files and chech that were succesfully closed
//       file.close();
//       if (file.is_open())
// 	throw 3;
//     }

//   // One more check, this time check if we store the same number of lines
//   // (just in case)
//   if (data.size() != date_and_time.size())
//     throw 4;    
//   // We're happy so far
//   // std::cout << "Done reading" << std::endl;
// }


void read_data (const std::vector< std::string >   &filenames,
		std::vector< std::vector<int> >    &date_and_time,
		std::vector< std::vector<double> > &data,
		bool day_number_column = true) // by default we assume that the file has a column for the day number
{
  /*
    We clear the containers that will store met
    data and date data from each set of files
  */
  data.clear();
  date_and_time.clear();
  /*
    In this loop, we will go through each set of 
    files, reading them, storing the entries in 
    temporal vectors and storing those vectors in
    the above containers
  */
  for (unsigned int j=0; j<filenames.size(); j++)
    {
      // Open files and check that were correctly opened
      std::ifstream file (filenames[j].c_str());
      if (!file.is_open())
	{
	  std::cout << "Error opening file: "
		    << filenames[j] << std::endl;
	  throw 2;
	}
      // Define some variables that will help us with the
      // reading process
      std::string line_file;
      std::string file_token;
      // We read each line in each file. If EOF is reach
      // in any of the files, the loop is finished
      while (std::getline(file,line_file)) 
	{
	  /*
	    there are two possible options (so far) for the file 
	    we are passing to this function. 1) the file is a met
	    file with date and time in the first and second column:
	    dd/mm/yyyy\tHH/MM/SS\tD1\tD2\t...
	    2) the file is a general file with only numbers in
	    columns:
	    D1\tD2\tD3\t.....
	    In the first case we look for the characters '\' or
	    ':' in the first line (we are assuming that all lines
	    follow the same pattern). If we find them then is the 
	    first case, if not, then is the second.
	  */
	  std::stringstream file_iss;
	  file_iss << line_file;
	  std::vector<int>    row_date_and_time;
	  std::vector<double> row_data;
	  if ((line_file.find('/') != std::string::npos) ||
	      (line_file.find(':') != std::string::npos))
	    {
	      std::stringstream   file_iss1;
	      int c = 0;
	  
	      while (std::getline(file_iss,file_token,'\t'))
		{
		  file_iss1 << file_token;
		  if (c==0)
		    while (std::getline(file_iss1,file_token,'/')) row_date_and_time.push_back(atoi(file_token.c_str()));
		  if (c==1)
		    while (std::getline(file_iss1,file_token,':')) row_date_and_time.push_back(atoi(file_token.c_str())); 
		  if (c>1)
		    while (std::getline(file_iss1,file_token,':'))
		      {
			if (day_number_column && c>2)
			  row_data.push_back(atof(file_token.c_str()));
			else if (!day_number_column)
			  row_data.push_back(atof(file_token.c_str()));
		      }
		  file_iss1.clear();
		  c++;
		}
	      file_iss.clear();
	      // Store temporal vectors to corresponding containers
	      date_and_time.push_back (row_date_and_time);
	      data.push_back          (row_data);
	    }
	  else
	    {
	      if (line_file.find('\t') != std::string::npos)
		while (std::getline(file_iss,file_token,'\t'))
		  row_data.push_back(atof(file_token.c_str()));
	      else
		while (std::getline(file_iss,file_token))
		  row_data.push_back(atof(file_token.c_str()));
	      
	      data.push_back(row_data);
	    }
	}
      // close files and check that were succesfully closed
      file.close();
      if (file.is_open())
	throw 3;
    }
}
//**************************************************************************************************
// Problem: The time seems to be shifted by 1 hour 
//
// Answer: saw  in StackOverflow: http://stackoverflow.com/questions/3660983/c-time-t-problem/3661129
//
// "mktime takes a struct tm giving a local time and returns the number of seconds since January 1, 1970,
// 0:00 UTC. Therefore, your GetDate(1970,1,1,0,0,0); call will return 0 if your local time zone is UTC,
// but may return other values for other time zones.
//
// Edit: For a UTC version of mktime or your GetDate, try the following (untested):
//    Call getenv to save the current value of the TZ environment variable (if any).
//    Call putenv to change the TZ environment variable to "UTC".
//    Call tzset to make your changes active.
//    Call mktime.
//    Restore the old value of TZ, then call _tzset again."
time_t mktimeUTC(struct tm* timeinfo)
{
  // *** enter in UTC mode
  char* oldTZ = getenv("TZ");
  char tzUTC[] = "TZ=UTC";
  putenv(tzUTC);
  tzset();
  // ***

  time_t ret = mktime ( timeinfo );

  // *** Restore previous TZ
  if(oldTZ == NULL)
    {
      char tz[] = "TZ=";
      putenv(tz);
    }
  else
    {
      char buff[255];
      sprintf(buff,"TZ=%s",oldTZ);
      putenv(buff);
    }
  tzset();
  // ***

  return ret;
}

//**************************************************************************************************
void date_to_seconds( const std::vector< std::vector<int> > &date_and_time,
		      std::vector< std::vector<int> > &date_in_seconds)
{
  time_t rawtime;// = time(0);
  time(&rawtime);
  struct tm *ref_date; // 01/07/2005   00:00:00
  
  // time_t rawtime = 0;
  // struct tm * ref_date; // 01/07/2005   00:00:00
  
  ref_date = localtime(&rawtime);  
  ref_date->tm_hour = 0;
  ref_date->tm_min  = 0;
  ref_date->tm_sec  = 0;
  ref_date->tm_mon  = 7-1;
  ref_date->tm_mday = 1;
  ref_date->tm_year = 2005-1970 + 70;
  time_t ref_seconds = mktimeUTC(ref_date);
  
  // std::cout << ref_date->tm_isdst << " " << ref_date->tm_yday << " " << ref_date->tm_wday << "\t"
  //        << ref_date->tm_mday  << "/" << ref_date->tm_mon  << "/" << ref_date->tm_year << "\t"
  //        << ref_date->tm_hour  << ":" << ref_date->tm_min  << ":" << ref_date->tm_sec  << "\t" 
  //        << ref_seconds  << std::endl;
  
  // Translate date and time entries to seconds since 01/07/2005
  for (unsigned int i=0; i<date_and_time.size(); i++)
    {
      time_t dummy_time = 0;
      // time(&dummy_time);
      struct tm *date;
      
      date = localtime(&dummy_time);
	
      date->tm_mday = date_and_time[i][0];
      date->tm_mon  = date_and_time[i][1]-1;
      date->tm_year = date_and_time[i][2]-1970+70;

      date->tm_hour = date_and_time[i][3];
      date->tm_min  = date_and_time[i][4];
      date->tm_sec  = date_and_time[i][5];
      date->tm_isdst = 1;
            
      std::vector<int> temp1;
      temp1.clear();
      temp1.push_back((int)difftime(mktimeUTC(date),ref_seconds));
      date_in_seconds.push_back(temp1);      
      // if (i<15)
      // 	std::cout << mktimeUTC(date) << "\t" << ref_seconds << "\t" << temp1[0] << std::endl;
      // if (i<15)
      // 	std::cout 
      // 	  << date->tm_mday << "/" << date->tm_mon << "/" << date->tm_year << "\t" 
      // 	  << date->tm_hour << ":" << date->tm_min << ":" << date->tm_sec << "\t" 
      // 	  << date_and_time[i][0] << "/" << date_and_time[i][1]-1 << "/" << date_and_time[i][2]-1900 << "\t" 
      // 	  << date_and_time[i][3] << ":" << date_and_time[i][4] << ":" << date_and_time[i][5] << "\t"
      // 	  << mktimeUTC(date) << "\t" << ref_seconds << "\t" << temp1[0] << std::endl;
    }
}
//**************************************************************************************************
void seconds_to_date( std::vector< std::vector<int> > &date_and_time,
		      const std::vector< std::vector<int> > &date_in_seconds)
{
  time_t rawtime = 0;
  struct tm * ref_date; // 01/07/2005   00:00:00
  ref_date = localtime(&rawtime);
  ref_date->tm_hour = 0;
  ref_date->tm_min  = 0;
  ref_date->tm_sec  = 0;
  ref_date->tm_mon  = 7-1;
  ref_date->tm_mday = 1;
  ref_date->tm_year = 2005-1970+70;

  time_t ref_seconds = mktimeUTC(ref_date);

  for (unsigned int i=0; i<date_in_seconds.size(); i++)
    {
      time_t time_entry = (time_t)date_in_seconds[i][0] + ref_seconds;
      struct tm * ptm;
      
      ptm = gmtime (&time_entry);
      std::vector<int> temp1;
      temp1.clear();
      temp1.push_back(ptm->tm_mday);
      temp1.push_back(ptm->tm_mon+1);
      temp1.push_back(ptm->tm_year+1900);
      temp1.push_back(ptm->tm_hour);
      temp1.push_back(ptm->tm_min);
      temp1.push_back(ptm->tm_sec);
      
      date_and_time.push_back(temp1);

      // if (i<15)
      // 	std::cout << ptm->tm_mday << "/" << ptm->tm_mon+1 << "/" << ptm->tm_year+1900 << "\t"
      // 		  << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << "\t"
      // 		  << date_in_seconds[i][0] << std::endl;
    }
}

//**************************************************************************************************
// linear interpolation
double interpolate_data(std::vector< std::pair<double,double> > table,
			const double x)
{
  const double INF = 1.e100;
  // Assumes that "table" is sorted by .first
  // Check if x is out of bound
  if (x > table.back().first)
    {
      std::cout << "Warning. In data_tools.h -> interpolate_data. value out of bound." << std::endl;
	return (INF);
    }
  if (x < table[0].first)
    {
      std::cout << "Warning. In data_tools.h -> interpolate_data. value out of bound." << std::endl;
      return (-INF);
    }
  std::vector<std::pair<double, double> >::iterator it, it2;
  // INFINITY is defined in math.h in the glibc implementation
  it = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -INF));
  // Corner case
  if (it == table.begin()) return (it->second);
  it2 = it;
  --it2;
  return (it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first));
}
//**************************************************************************************************
void make_scatter_plots(const std::string script_filename,
			const std::string data_output_filename,
			const std::vector< std::string > graph_filenames,
			const std::vector< double > &xdata,
			const std::vector< double > &ydata,
			const std::vector< std::vector<int> > &date_and_time,
			std::string main_variable_name,
			std::string main_variable_units,
			std::string optional_variable_name = "",
			std::string optional_variable_units = "",
			double optional_variable_value = 0.)
{
  if(xdata.size()!=ydata.size())
    throw 4;
  /*
    Generate GNUplot scripts to plot the actual
    scatter graphs.
    First we calculate the coefficients for the 
    least squares line adjusment. The idea is to
    include in the scatter graph a line of the form
    Y = Ax + B
    The following code is to calculate A and B.
  */
  double Xmax  = 0;
  double Xmin  = 10000;
  double Ymax  = 0;
  double Ymin  = 10000;
  double Xmean = 0;
  double Ymean = 0;
  double XYmean = 0;
  double X2mean = 0;
  double Y2mean = 0;
    
  double MSE   = 0;
  double RMSD  = 0;
  double NRMSD = 0;
  
  for (unsigned int i=0; i<xdata.size(); i++)
    {
      double x = xdata[i];
      double y = ydata[i];
      
      if (x > Xmax)
  	Xmax = x;
      if (x < Xmin)
  	Xmin = x;
      if (y > Ymax)
  	Ymax = y;
      if (y < Ymin)
  	Ymin = y;
      
      Xmean  += x;
      Ymean  += y;
      XYmean += x*y;
      X2mean += x*x;
      Y2mean += y*y;
	
      MSE += pow(x-y,2); 
    }
  
  int n = xdata.size();

  Xmean = Xmean/n;
  Ymean = Ymean/n;
  XYmean = XYmean/n;
  X2mean = X2mean/n;
  Y2mean = Y2mean/n;
  
  MSE   = MSE/n;
  RMSD  = sqrt(MSE);
  NRMSD = RMSD/(Xmax-Xmin);

  double alpha  = 0;
  double beta   = 0;
  double gamma  = 0;
  double beta_1 = 0;
  double beta_2 = 0;
  double SS_res = 0;
  double SS_tot = 0;
  
  // Fitting the regression line  Y = alpha + beta*X  AND  Y = gamma*X
  for (unsigned int i=0; i<xdata.size(); i++)
    {
      double x = xdata[i];
      double y = ydata[i];
	
      beta_1 += (x-Xmean)*(y-Ymean);
      beta_2 += pow(x-Xmean,2);
    }
  beta  = beta_1/beta_2;
  alpha = Ymean - beta*Xmean;
  gamma = XYmean/X2mean;

  double Fmax = 0;
  double Fmin = 0;
  for (unsigned int i=0; i<xdata.size(); i++)
    {
      double x = xdata[i];
      double y = ydata[i];
      double f = alpha + beta*x;

      if (f > Fmax)
  	Fmax = f;
      if (f < Fmin)
  	Fmin = f;

      SS_tot += pow(y-Ymean,2); //the total sum of squares (proportional to the sample variance)      
      SS_res += pow(y-f,2);     //the sum of squares of residuals, also called the residual sum of squares.
    }
  double R2 = 0;
  R2 = 1-SS_res/SS_tot;
  /*
    Generate files with data to draw scatter plots    
    date     time     file1_entry     file2_entry     difference     normalized difference     day     month     year
  */
  std::vector< std::vector<double> > processed_data;
  processed_data.clear();
  for (unsigned int i=0; i<xdata.size(); i++)
    {
      std::vector<double> temp;
      temp.clear();
      temp.push_back(xdata[i]);
      temp.push_back(ydata[i]);
      temp.push_back(ydata[i] - (alpha + beta*xdata[i]));
      temp.push_back(pow(ydata[i] - (alpha + beta*xdata[i]),2));
      temp.push_back(pow(ydata[i]-Ymean,2));
      temp.push_back((double) date_and_time[i][0]);
      temp.push_back((double) date_and_time[i][1]);
      temp.push_back((double) date_and_time[i][2]);

      processed_data.push_back(temp);
    }
  
  std::ofstream data_output_file (data_output_filename.c_str());
  if (!data_output_file.is_open())
    throw 2;
  print_data (0,
	      date_and_time,
	      processed_data,
	      /*std::cout*/data_output_file);
  data_output_file.close();
  if (data_output_file.is_open())
    throw 3;

  // Now, let us write the actual script for gnuplot.
  // We make use of AWK program. The idea is to plot
  // with a different color the data for corresponding
  // to different months

  std::ofstream outFile (script_filename.c_str());
  if (!outFile.is_open())
    throw 2;
  
  outFile << "#*********************************************************"
  	  << std::endl
  	  << " \t #name of output file: "  << data_output_filename << std::endl
  	  << std::endl
  	  << "\t #Max X value  : " << Xmax  << std::endl
  	  << "\t #Min X value  : " << Xmin  << std::endl
  	  << "\t #Mean X value : " << Xmean << std::endl
  	  << "\t #Max Y value    : " << Ymax  << std::endl
  	  << "\t #Min Y value    : " << Ymin  << std::endl
  	  << "\t #Mean Y value   : " << Ymean << std::endl
  	  << "\t #MSE   : " << MSE   << std::endl
  	  << "\t #RMSD  : " << RMSD  << std::endl
  	  << "\t #NRMSD : " << NRMSD << std::endl
  	  << "\t #R2    : " << R2    << std::endl
  	  << std::endl
  	  << "\t #Equation: \t" << "Xvalue = (" << alpha << ") + (" << beta << ")*Yvalue" 
  	  << std::endl
  	  << "#*********************************************************"
  	  << std::endl;
   
  outFile << "a =" << alpha << std::endl
  	  << "b =" << beta  << std::endl
  	  << "c =" << gamma << std::endl
  	  << "set title " << "\"" << main_variable_name << " -  Experimental vs Predicted";
  
  if (optional_variable_value != 0.)
    outFile << "\\n" << optional_variable_name << " " << optional_variable_value << " " << optional_variable_units;
  outFile << "\"" << std::endl
  	  << "set xlabel \"Experimental data " << main_variable_units << "\"" << std::endl
  	  << "set ylabel \"Predicted data "   << main_variable_units << "\"" << std::endl
  	  << "set size square"   << std::endl;
  
  double range_max = 0;
  double range_min = 0;
  double tics = 0;
  double Max = 0;
  double Min = 0;

  if (Xmax>Ymax)
    Max=Xmax;
  else
    Max=Ymax;

  if (Xmin<Ymin)
    Min=Xmin;
  else
    Min=Ymin;

  if ((Max-Min)<10)
    {
      range_max = ceil  (Max);
      range_min = floor (Min);
      tics = 0.5;
    }
  else if ((Max-Min)<100)
    {
      range_max = 10.*ceil  (Max/10);
      range_min = 10.*floor (Min/10);
      tics = 2;
    }
  else if ((Max-Min)<500)
    {
      range_max = 10.*ceil  (Max/10);
      range_min = 10.*floor (Min/10);
      tics = 50;
    }
  else 
    {
      range_max = 100.*ceil  (Max/100);
      range_min = 100.*floor (Min/100);
      tics = 100;
    }
  /*
    Scatter plot
  */  
  outFile << "set terminal pngcairo dashed size 1700,1000 font \"Arial, 18\"" << std::endl
  	  << "set output \"" << graph_filenames[0] << "\"" << std::endl;
  outFile << "set xrange [" << range_min << ":" << range_max << "]" << std::endl
  	  << "set yrange [" << range_min << ":" << range_max << "]" << std::endl
  	  << "set xtics " << tics << std::endl
  	  << "set ytics " << tics << std::endl
  	  << "set grid"           << std::endl
  	  << "set key vert at " << range_max + tics/5 << "," << range_max << " left Left" << std::endl
  	  << "set label 1 \"********************\\n RMSD = " << RMSD << main_variable_units << "\\n NRMSD = " << NRMSD << "\\n R2    = " << R2 
	  << " \\n********************\" at " << range_max + tics/5 << "," << 0.5*(range_max+range_min) - tics << std::endl
	  << "set lmargin 0" << std::endl
  	  << std::endl;
  outFile << "plot\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"9.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt  9 pt  9 t \"September\",\\" << std::endl
  	  << "\"< awk '{if(($9 == \\\"10.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 10 pt 10 t \"October  \",\\" << std::endl
  	  << "\"< awk '{if(($9 == \\\"11.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 11 pt 11 t \"November \",\\" << std::endl
  	  << "\"< awk '{if(($9 == \\\"12.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 12 pt 12 t \"December \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"1.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 13 pt 13 t \"January  \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"2.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 14 pt 14 t \"February \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"3.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 15 pt 15 t \"March    \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"4.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 16 pt 16 t \"April    \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"5.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 17 pt 17 t \"May      \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"6.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 18 pt 18 t \"June     \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"7.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 19 pt 19 t \"July     \",\\" << std::endl
	  << "\"< awk '{if(($9 ==  \\\"8.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 20 pt 20 t \"August   \",\\" << std::endl
      	  << "x"       << " lt 2 lc 4 lw 1 t \"Ideal\",\\" << std::endl
  	  << "a + b*x" << " lt 1 lc 7 lw 2 t \"Y = " << alpha << "\\n    + " << beta << "*X\"" 
  	  << std::endl
  	  << std::endl;
  
  /*
    Residual Plot
  */
  if (Fmax>Ymax)
    Max=Fmax;
  else
    Max=Ymax;

  if (Fmin<Ymin)
    Min=Fmin;
  else
    Min=Ymin;
  
  if ((Max-Min)<10)
    {
      range_max = ceil  (Max);
      range_min = floor (Min);
      tics = 0.5;
    }
  else if ((Max-Min)<100)
    {
      range_max = 10.*ceil  (Max/10);
      range_min = 10.*floor (Min/10);
      tics = 2;
    }
  else if ((Max-Min)<500)
    {
      range_max = 100.*ceil  (Max/100);
      range_min = 100.*floor (Min/100);
      tics = 50;
    }
  else 
    {
      range_max = 100.*ceil  (Max/100);
      range_min = 100.*floor (Min/100);
      tics = 100;
    }

  outFile << "set terminal pngcairo size 1300,1000 font \"Arial, 18\"" << std::endl
  	  << "set output \"" << graph_filenames[1] << "\"" << std::endl;

  outFile << "set title " << "\"" << main_variable_name << " - Residuals";
  if (optional_variable_value != 0.)
    outFile << "\\n" << optional_variable_name << " " << optional_variable_value << " " << optional_variable_units;
  outFile << "\"" << std::endl;
  
  outFile << "set xrange [" << range_min << ":" << range_max << "]" << std::endl
  	  << "set yrange [" << range_min << ":" << range_max << "]" << std::endl
	  << "unset label 1" << std::endl
  	  << "set key vert at " << range_max + tics/5 << "," << range_max << " left Left" << std::endl
  	  << "set grid"           << std::endl
	  << "set xlabel \"Predicted data " << main_variable_units << "\"" << std::endl
  	  << "set ylabel \"Residuals "       << main_variable_units << "\"" << std::endl
	  
	  << "set lmargin 0" << std::endl
	  << std::endl; 
  
  outFile << "plot\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"9.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt  9 pt  9 t \"September\",\\" << std::endl
  	  << "\"< awk '{if(($9 == \\\"10.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 10 pt 10 t \"October  \",\\" << std::endl
  	  << "\"< awk '{if(($9 == \\\"11.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 11 pt 11 t \"November \",\\" << std::endl
  	  << "\"< awk '{if(($9 == \\\"12.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 12 pt 12 t \"December \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"1.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 13 pt 13 t \"January  \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"2.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 14 pt 14 t \"February \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"3.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 15 pt 15 t \"March    \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"4.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 16 pt 16 t \"April    \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"5.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 17 pt 17 t \"May      \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"6.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 18 pt 18 t \"June     \",\\" << std::endl
  	  << "\"< awk '{if(($9 ==  \\\"7.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 19 pt 19 t \"July     \",\\" << std::endl
	  << "\"< awk '{if(($9 ==  \\\"8.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 20 pt 20 t \"August   \"\\" << std::endl	
  	  << std::endl
  	  << std::endl;
  /*
    Plot time series
  */
  outFile << "set terminal pngcairo size 1700,1000 font \"Arial, 18\"" << std::endl
  	  << "set output \"" << graph_filenames[2] << "\"" << std::endl;
  
  outFile << "set title " << "\"" << main_variable_name << " - Time series";
  if (optional_variable_value != 0.)
    outFile << "\\n" << optional_variable_name << " " << optional_variable_value << " " << optional_variable_units;
  outFile << "\"" << std::endl;

  outFile << "set xdata time" << std::endl
	  << "set timefmt \"%d/%m/%Y\\t%H:%M:%S\"" << std::endl
	  << "set format x \"%d/%m/%y\"" << std::endl
	  << "set xtics auto" << std::endl
	  << "set size nosquare"   << std::endl;
  
  outFile << "set autoscale" << std::endl
	  << "set yrange [" << range_min << ":" << range_max << "]" << std::endl
	  << "set key default" << std::endl
  	  << "set grid"           << std::endl
	  << "set xlabel \"Time\"" << std::endl
  	  << "set ylabel \"" << main_variable_units << "\"" << std::endl
	  << "set lmargin" << std::endl
	  << std::endl;
  
  outFile << "plot\\" << std::endl
	  << "\"" << data_output_filename << "\" u 1:4 w lp lt 25 t \"Predicted   \",\\" << std::endl
	  << "\"" << data_output_filename << "\" u 1:3 w lp lt 1  t \"Experimental\"\\" << std::endl;

  outFile.close();
  if (outFile.is_open())
    throw 3;


}

//*********************************************************************************
//*********************************************************************************
//*********************************************************************************

class Names
{
public:
  Names (const int preheating_step_,
	 const std::string met_data_type_);
  std::vector<double> soil_depths;
  std::vector<double> road_depths;
  std::vector<std::string> bh_temperatures;
  std::vector<std::string> trl_met_data;
  std::vector<std::string> met_office_daily_averages;
private:
  const int preheating_step;
  const std::string met_data_type;
};

Names::Names (const int preheating_step_,
	      const std::string met_data_type_):
  preheating_step (preheating_step_),
  met_data_type   (met_data_type_)
{
  soil_depths.push_back(  0.025 );
  soil_depths.push_back(  0.125 );
  soil_depths.push_back(  0.825 );
  soil_depths.push_back(  0.875 );
  soil_depths.push_back(  1.025 );
  soil_depths.push_back(  1.175 );
  soil_depths.push_back(  1.375 );
  soil_depths.push_back(  1.875 );
  soil_depths.push_back(  3.875 );
  soil_depths.push_back(  7.875 );
  soil_depths.push_back( 12.875 );

  road_depths.push_back(  0.000); // 0
  road_depths.push_back(  0.010); // 0
  road_depths.push_back(  0.025); // 1
  road_depths.push_back(  0.050); // 2
  road_depths.push_back(  0.075); // 3
  road_depths.push_back(  0.100); // 4

  road_depths.push_back(  0.120); // 5
  road_depths.push_back(  0.1325); // 6
  road_depths.push_back(  0.160); // 7
  road_depths.push_back(  0.180); // 8

  road_depths.push_back(  0.200); // 9
  road_depths.push_back(  0.250); // 10
  road_depths.push_back(  0.300); // 11
  road_depths.push_back(  0.350); // 12
  road_depths.push_back(  0.400); // 13
  road_depths.push_back(  0.450); // 14
  road_depths.push_back(  0.500); // 15
  road_depths.push_back(  0.550); // 16
  road_depths.push_back(  0.600); // 17
  road_depths.push_back(  0.650); // 18
  road_depths.push_back(  0.700); // 19
  road_depths.push_back(  0.750); // 20
  road_depths.push_back(  0.800); // 21
  road_depths.push_back(  0.8475); // 22

  road_depths.push_back(  0.875); // 23
  road_depths.push_back(  0.925); // 24
  road_depths.push_back(  0.975); // 25
  road_depths.push_back(  1.025); // 26
  road_depths.push_back(  1.175); // 27
  road_depths.push_back(  1.375); // 28
  road_depths.push_back(  1.875); // 29
  road_depths.push_back(  2.875); // 30
  road_depths.push_back(  3.875); // 31
  road_depths.push_back(  7.875); // 32
  road_depths.push_back( 12.875); // 33
  
  bh_temperatures.push_back ("./input/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_09.txt");
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_10.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_11.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_12.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_01.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_02.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_03.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_04.txt"); */
  /* bh_temperatures.push_back ("./input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_05.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_06.txt"); */
  /* bh_temperatures.push_back ("/home/zerpiko/Dropbox/PhD/input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_07.txt"); */
  /* bh_temperatures.push_back ("./input/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_08.txt"); */
    
  
  if (met_data_type=="met_office_data")
    {
      /*
	Metereological data from Met Office
      */
      if (preheating_step==1)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_09.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_10_using_2006_data.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_11.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_12.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_01.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_02.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_03.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_04.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_05.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_06.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_07.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_08.txt");
	}
      else if (preheating_step==2)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_09.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_10_using_2006_data.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_11.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2005_12.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_01.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_02.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_03.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_04.txt");
	}
      else if (preheating_step==3)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_05.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_06.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_07.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_08_up_to_day_22.txt");
	}
      else if (preheating_step==4)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/hourly_average_met_office_data_from_23_08_2005_to_14_11_2005.txt");
	}
      else
	{
	  std::cout << "Error in Names class. Meteorological data from Met Office." 
		    <<  "Preheating step not implemented." << std::endl;
	  throw 1;
	}
    }
  else if (met_data_type=="trl_met_data")
    {
      /*
	Metereological data from TRL
      */
      if (preheating_step==1)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_09.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_10.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_11.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_12.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_01.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_02.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_03.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_04.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_05.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_06.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_07.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_08.txt");
	}
      else if (preheating_step==2)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_09.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_10.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_11.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_05_12.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_01.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_02.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_03.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_04.txt");
	}
      else if (preheating_step==3)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_05.txt"); 
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_06.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_06_07.txt");
	  trl_met_data.push_back ("/home/c1045890/input/met_data/met_office/met_office_2006_08_up_to_day_22.txt");
	}
      else if (preheating_step==4)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_from_23_08_2005_to_14_11_2005.txt");
	}
      else if (preheating_step==5)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_from_15_11_2005_to_20_02_2006.txt");
	}
      else if (preheating_step==6)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_from_21_02_2006_to_26_04_2006.txt");
	}
      else if (preheating_step==7)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_from_27_04_2006_to_31_10_2006.txt");
	}
      else if (preheating_step==8)
	{
	  trl_met_data.push_back ("/home/c1045890/input/met_data/trl/hourly_average_met_data_from_01_11_2006_to_28_02_2007.txt");
	}
      else
	{
	  std::cout << "Error in Names class. Meteorological data from TRL." 
		    <<  "Preheating step not implemented." << std::endl;
	  throw 1;
	}
    }
  else
    {
      std::cout << "Error in Names class.\n" 
		<< "Met data type not implemented.\n"
		<< "Options are:\n"
		<< "\t --\"trl_met_data\"\n"
		<< "\t --\"met_office_data\"\n" << std::endl;
      throw 1;
    }
  
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1985.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1986.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1987.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1988.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1989.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1990.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1991.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1992.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1993.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1994.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1995.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1996.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1997.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1998.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_1999.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2000.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2001.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2002.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2003.dat");
  met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2004.dat");
  // met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2005.dat");
  // met_office_daily_averages.push_back("/home/c1045890/input/met_data/met_office/daily_averages_solar_air_temperature_2006.dat");
}

//*********************************************************************************
//*********************************************************************************
//*********************************************************************************

void read_met_data (  std::vector< std::vector<int> >    &date_and_time,
		      std::vector< std::vector<double> > &met_data,
		      int time_step,
		      const int preheating_step,
		      const std::string met_data_type)
{
  date_and_time.clear();
  met_data.clear();
  
  bool day_number = false;  
  
  std::string met_data_filename = "met_data_ph_step_";
  std::stringstream ph_step;
  ph_step << preheating_step;
  met_data_filename += ph_step.str() + ".txt";
    
  /*
    Check if we have a file with the previous name. If so
    read from this file the met data. We are assuming that 
    the file contains the corresponding met data for the 
    preheating step of interest. There is not way to check
    if this actually true. I'm just guessing that there
    should be any other file there that happens to have
    the same name and doesn't contain the desired met data.
  */
  std::ifstream file (met_data_filename.c_str());
  if (file.good())
    {
      std::vector<std::string> met_data_file;
      met_data_file.push_back(met_data_filename);
      
      read_data (met_data_file,
		 date_and_time,
		 met_data,
		 day_number);
    }
  else
    {
      Names names (preheating_step,
		   met_data_type);
      read_data (names.trl_met_data,
		 date_and_time,
		 met_data,
		 day_number);
      if (time_step<3600)
	{
	  std::vector< std::vector<int> > date_and_time_in_seconds;
      
	  std::vector< std::vector<int> >    temp_date_and_time;
	  std::vector< std::vector<int> >    temp_date_and_time_in_seconds;
	  std::vector< std::vector<double> > temp_met_data;
      
	  date_to_seconds(date_and_time,
			  date_and_time_in_seconds);
	  if (date_and_time_in_seconds.size() != date_and_time.size())
	    throw 3;
	
	  // generate the output times we want
	  int seconds = date_and_time_in_seconds[0][0];
	  while (seconds <= date_and_time_in_seconds.back()[0])
	    {
	      std::vector<int> temp_date_and_time_in_seconds_row;
	      temp_date_and_time_in_seconds_row.push_back (seconds);
	      seconds += time_step;
	      temp_date_and_time_in_seconds.push_back(temp_date_and_time_in_seconds_row);
	    }
	  // generate the corresponding dates
	  seconds_to_date( temp_date_and_time,
			   temp_date_and_time_in_seconds);
	
	  std::vector<std::vector<std::pair<double,double> > >tables;
	  for (unsigned int j=0; j<met_data[0].size(); j++)
	    {
	      std::vector<std::pair<double,double> > table;
	      for (unsigned int i=0; i<met_data.size(); i++)
		table.push_back(std::make_pair((double)date_and_time_in_seconds[i][0],met_data[i][j]));
	      tables.push_back(table);
	    }
	  // interpolate data according to the previous
	  // generated time entries
	  for (unsigned int i=0; i<temp_date_and_time_in_seconds.size(); i++)
	    {
	      std::vector<double> temp1;
	      temp1.clear();
	      for (unsigned int j=0; j<met_data[0].size(); j++)
		temp1.push_back( interpolate_data(tables[j],(double)temp_date_and_time_in_seconds[i][0]));
	      temp_met_data.push_back(temp1);
	    }

	  met_data = temp_met_data;
	  date_and_time = temp_date_and_time;
	}
      /*
	If we need to consider several years we do it like this
      */
      if (preheating_step==1)
	{
	  unsigned int one_year_size = date_and_time.size();
	  for  (unsigned int i=0; i<7; i++)
	    for (unsigned int j=0; j<one_year_size; j++)
	  {
	    std::vector<int> date;
	    date=date_and_time[j];
	    date[2]= date[2]+i+1;
	    
	    date_and_time.push_back(date);
	    met_data.push_back(met_data[j]);
	  }
	}
      
      /*
	Print met data file
      */
      {
        std::ofstream file (met_data_filename.c_str());
        if (!file.is_open())
          throw 2;
	
	print_data(0,
		   date_and_time,
		   met_data,
		   file);
	
        file.close();
        if (file.is_open())
          throw 3;
      }
    }
}
