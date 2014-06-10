//-*-C++-*-

#ifndef TIME_DOMAIN_H
#define TIME_DOMAIN_H

/*!
 *  \author Peter Staar
 */
class time_domain 
{
public:

  const static int RULE      = 1;
  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  typedef MATH_ALGORITHMS::interval_dmn_1D_type dmn_specifications_type;

public:

  //const static int DIMENSION = 1;
  typedef time_domain parameter_type; // --> used in the interpolation!

  static int&                 get_size();

  static std::string          get_name();
  
  static std::vector<double>& get_elements();

  static double get_beta();

  static double get_volume();

  template<class stream_type>
  static void to_JSN(stream_type& ss);

  template<IO::FORMAT DATA_FORMAT>
  static void read(IO::reader<DATA_FORMAT>& reader);

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& writer);

  template<typename parameters_t>
  static void initialize(parameters_t& parameters);

  static void initialize_integration_domain(int                        level,
					    std::vector<scalar_type>&  weights,
					    std::vector<element_type>& elements);

private:

  static int    time_slices;
  static double beta;
};

int    time_domain::time_slices = -1;
double time_domain::beta        = 0;

int& time_domain::get_size()
{
  static int size = -1;
  return size;
}

std::string time_domain::get_name()
{
  static std::string name = "time-domain";
  return name;
}

std::vector<double>& time_domain::get_elements()
{
  static std::vector<double> elements;
  return elements;
}

double time_domain::get_beta(){
  return (get_elements().back()+1.e-10);
}

double time_domain::get_volume(){
  return (get_elements().back()+1.e-10);
}

template<IO::FORMAT DATA_FORMAT>
void time_domain::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(get_name());

  writer.execute("elements", get_elements());

  writer.close_group();
}

template<typename parameters_t>
void time_domain::initialize(parameters_t& parameters)
{
  //time_slices = parameters.get_number_of_positive_times();
  time_slices = parameters.get_sp_time_intervals();
  beta        = parameters.get_beta();

  get_size() = 2*(parameters.get_sp_time_intervals()+1);

  get_elements().resize(get_size());

  for(int i=0; i<get_size()/2; i++){
    get_elements()[i+get_size()/2] = double(i)/(double(get_size())/2.-1.)*beta;
    get_elements()[i             ] = -beta + double(i)/(double(get_size())/2.-1.)*beta;
  }
  
  get_elements()[0]              += 1.e-10;
  get_elements()[get_size()/2-1] -= 1.e-10;
  get_elements()[get_size()/2]   += 1.e-10;
  get_elements()[get_size()/2-1] -= 1.e-10;
}

template<class stream_type>
void time_domain::to_JSN(stream_type& ss)
{
  ss << "\"time_domain\" : [\n";

  for(int i=0; i<get_size(); i++)
    if(i == get_size()-1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";

  ss << "]\n";
}

/*
void time_domain::initialize_integration_domain(int                        level,
						std::vector<scalar_type>&  weights,
						std::vector<element_type>& elements)
{
  std::vector<std::pair<double,double> > vecs(0);

  {
    int point_num = 0;
    MATH_ALGORITHMS::GAUSSIAN_QUADRATURE::get_size(RULE, DIMENSION, point_num);

    scalar_type* w = new scalar_type[point_num];
    scalar_type* x = new scalar_type[point_num*DIMENSION];

    MATH_ALGORITHMS::GAUSSIAN_QUADRATURE::get_points(RULE, DIMENSION, point_num, w, x);

    for(int j=0; j<point_num; j++)
      {
	std::pair<double,double> p(x[j], w[j]);
	vecs.push_back(p);
      }

    sort(vecs.begin(), vecs.end());

    for(int j=0; j<vecs.size()-1; j++){
      if(abs(vecs[j+1].first-vecs[j].first)<1.e-6){
	vecs[j].second += vecs[j+1].second;
	
	vecs.erase(vecs.begin()+j+1);
      }
    }

    delete [] w;
    delete [] x;
  
  }

  {
    int N = std::pow(2., level);
    
    weights .resize(0);
    elements.resize(0);
    
    for(int i=0; i<get_size()/2-2; i++){

      element_type min = get_elements()[i];
      element_type max = get_elements()[i+1];
      
      double delta = double(max-min)/double(N);
      
      for(int j=0; j<N; j++){
	  
	for(int l=0; l<vecs.size(); l++){
	  
	  elements.push_back(min + j*delta*vecs[l].first);
	  
	  weights .push_back(vecs.size()*vecs[l].second*get_volume());
	}
      }
    }
    
    for(int i=get_size()/2; i<get_size()-1; i++){

      element_type min = get_elements()[i];
      element_type max = get_elements()[i+1];
      
      double delta = double(max-min)/double(N);
	
      for(int j=0; j<N; j++){
	
	for(int l=0; l<vecs.size(); l++){
	  
	  elements.push_back(min + j*delta*vecs[l].first);
	  
	  weights .push_back(vecs.size()*vecs[l].second*get_volume());
	}
      }
    }
  }
}
*/

void time_domain::initialize_integration_domain(int                        level,
						std::vector<scalar_type>&  weights,
						std::vector<element_type>& elements)
{
  int N = std::pow(2., level);
    
  weights .resize(0);
  elements.resize(0);
    
//   {
//       elements.push_back(get_elements()[0]);      
//       weights .push_back(get_volume()/2.);
//   }

  for(int i=0; i<get_size()/2-1; i++){
      
    element_type min = get_elements()[i];
    element_type max = get_elements()[i+1];
    
    double delta = double(max-min)/double(N);
    
    for(int j=0; j<N; j++){
      
      elements.push_back(min + j*delta);
      
      weights .push_back(get_volume());
    }
  }

//   {
//       elements.push_back(get_elements()[get_size()/2-1]);      
//       weights .push_back(get_volume()/2.);
//   }

//   {
//       elements.push_back(get_elements()[get_size()/2]);      
//       weights .push_back(get_volume()/2.);
//   }
  
  for(int i=get_size()/2; i<get_size()-1; i++){
    
    element_type min = get_elements()[i];
    element_type max = get_elements()[i+1];
    
    double delta = double(max-min)/double(N);
    
    for(int j=0; j<N; j++){
      
      elements.push_back(min + j*delta);
      
      weights .push_back(get_volume());
    }

//     {
// 	elements.push_back(get_elements()[get_size()-1]);      
// 	weights .push_back(get_volume()/2.);
//     }
  }  
}

#endif




