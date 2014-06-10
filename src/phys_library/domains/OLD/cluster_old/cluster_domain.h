//-*-C++-*-

#ifndef CLUSTER_DOMAIN_H
#define CLUSTER_DOMAIN_H

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
struct cluster_specifications_type
{};

template<typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications_type<scalar_type, CLUSTER, MOMENTUM_SPACE, S>
{
  typedef MATH_ALGORITHMS::domain_specifications<scalar_type, std::vector<scalar_type>, 
						 MATH_ALGORITHMS::DISCRETE, MATH_ALGORITHMS::KRONECKER_DELTA, 
						 MATH_ALGORITHMS::PERIODIC, MATH_ALGORITHMS::EQUIDISTANT    > dmn_specifications_type;
};

template<typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications_type<scalar_type, CLUSTER, REAL_SPACE, S>
{
  typedef MATH_ALGORITHMS::domain_specifications<scalar_type, std::vector<scalar_type>, 
						 MATH_ALGORITHMS::EXPANSION, MATH_ALGORITHMS::HARMONICS, 
						 MATH_ALGORITHMS::PERIODIC , MATH_ALGORITHMS::EQUIDISTANT> dmn_specifications_type;
};

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class cluster_domain
{
public:

  static int DIMENSION;// = -1;

  const static CLUSTER_NAMES          NAME           = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE          SHAPE          = S;

  typedef typename cluster_specifications_type<scalar_type, N, R, S>::dmn_specifications_type dmn_specifications_type;

  //typedef typename dmn_specifications_type::scalar_type  scalar_type;
  //typedef typename dmn_specifications_type::element_type element_type;

  typedef std::vector<scalar_type> element_type;

public:
  
  static bool& is_initialized();
  
  static int& get_size();

  static int*& get_dimensions();

  static scalar_type*& get_basis();
  static scalar_type*& get_super_basis();

  static std::string& get_name();

  static std::vector<element_type>& get_elements();

  static scalar_type volume();
  static int         origin_index();
  
  static int add     (int i, int j);
  static int subtract(int i, int j);

  static LIN_ALG::matrix<int, LIN_ALG::CPU>& get_add_matrix();
  static LIN_ALG::matrix<int, LIN_ALG::CPU>& get_subtract_matrix();

  static void reset();

  template<typename ss_type>
  static void print(ss_type& ss);

  template<class stream_type>
  static void to_JSON(stream_type& ss);

};

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, N, R, S>::DIMENSION = -1;

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
bool& cluster_domain<scalar_type, N, R, S>::is_initialized()
{
  static bool initialized = false;
  return initialized;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int& cluster_domain<scalar_type, N, R, S>::get_size()
{
  static int size = 0;
  return size;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int*& cluster_domain<scalar_type, N, R, S>::get_dimensions()
{
  static int* dimensions = NULL;
  return dimensions;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& cluster_domain<scalar_type, N, R, S>::get_basis()
{
  static scalar_type* basis = NULL;
  return basis;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& cluster_domain<scalar_type, N, R, S>::get_super_basis()
{
  static scalar_type* basis = NULL;
  return basis;
}
  
template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::string& cluster_domain<scalar_type, N, R, S>::get_name()
{
  static std::string name = to_str(NAME)+" "+to_str(REPRESENTATION)+" "+to_str(SHAPE)+" (DIMENSION : " + to_str(DIMENSION) + ")";
  return name;
}
  
template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type> >& cluster_domain<scalar_type, N, R, S>::get_elements()
{
  static std::vector<element_type> elements(get_size());
  return elements;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type cluster_domain<scalar_type, N, R, S>::volume()
{
  static scalar_type volume = 0;
  return volume;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, N, R, S>::origin_index()
{
  static int index = cluster_operations::origin_index(get_elements());
  assert(VECTOR_OPERATIONS::L2_NORM(get_elements()[index])<1.e-6);
  return index;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, N, R, S>::add(int i, int j)
{
  static LIN_ALG::matrix<int, LIN_ALG::CPU>& A = get_add_matrix();
  return A(i,j);
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
LIN_ALG::matrix<int, LIN_ALG::CPU>& cluster_domain<scalar_type, N, R, S>::get_add_matrix()
{
  static LIN_ALG::matrix<int, LIN_ALG::CPU> A("add", get_size());
  return A;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, N, R, S>::subtract(int i, int j)
{
  static LIN_ALG::matrix<int, LIN_ALG::CPU>& A = get_subtract_matrix();
  return A(i,j);
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
LIN_ALG::matrix<int, LIN_ALG::CPU>& cluster_domain<scalar_type, N, R, S>::get_subtract_matrix()
{
  static LIN_ALG::matrix<int, LIN_ALG::CPU> A("subtract", get_size());
  return A;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
void cluster_domain<scalar_type, N, R, S>::reset()
{
  get_size() = 0;
  get_name() = "";
  
  get_elements().resize(0);
  
  is_initialized() = false;
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
template<typename ss_type>
void cluster_domain<scalar_type, N, R, S>::print(ss_type& ss)
{
  ss << get_name() << "\n";
  ss << "\t origin-index : " << origin_index() << "\n";
 
  ss << "\t basis : \n";
  for(int d0=0; d0<DIMENSION; d0++){
    ss << "\t";
    for(int d1=0; d1<DIMENSION; d1++)
      ss << get_basis()[d0+d1*DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t super-basis : \n";
  for(int d0=0; d0<DIMENSION; d0++){
    ss << "\t";
    for(int d1=0; d1<DIMENSION; d1++)
      ss << get_super_basis()[d0+d1*DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  for(int l=0; l<get_size(); l++){
    ss << "\t" << l << "\t";
    VECTOR_OPERATIONS::PRINT(get_elements()[l]);
    ss << "\n";
  }    

  get_add_matrix().print();

  get_subtract_matrix().print();
}

template<typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
template<typename ss_type>
void cluster_domain<scalar_type, N, R, S>::to_JSON(ss_type& ss)
{

}

#endif
