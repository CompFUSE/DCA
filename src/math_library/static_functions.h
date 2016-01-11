//-*-C++-*-

double sgn(double a)
{ 
  if(a>=0)
    return 1.;
  else
    return -1.;
}

template<typename scalartype>
inline scalartype abs_value(scalartype a)
{ 
  return std::fabs(a); 
}

template<typename scalartype>
inline scalartype abs_value(std::complex<scalartype> a)
{ 
  return abs(a); 
}

template<typename scalartype>
inline scalartype conj_value(scalartype a)
{ 
  return a; 
}

template<typename scalartype>
inline std::complex<scalartype> conj_value(std::complex<scalartype> a)
{ 
  return conj(a); 
}

template<typename scalartype>
scalartype square(scalartype a)
{ 
  return a*a; 
}

template<typename scalartype>
double cot(scalartype z) 
{ 
  return 1.0 / tan(z); 
} 

template<typename scalartype>
double sec(scalartype z) 
{ 
  return 1.0 / cos(z); 
} 

template<typename scalartype>
double csc(scalartype z) 
{ 
  return 1.0 / sin(z); 
} 

inline double erf_inverse_help(double t)
{
  // Abramowitz and Stegun formula 26.2.23.
  // The absolute value of the error should be less than 4.5 e-4.
  static double c[3] = {2.515517, 0.802853, 0.010328};
  static double d[3] = {1.432788, 0.189269, 0.001308};
  
  return t - ((c[2]*t + c[1])*t + c[0]) / 
    (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

inline double erf_inverse(double p)
{
  assert(p>1.e-16 and p<1-1.e-16);

  // See article above for explanation of this section.
  if(p < 0.5)
    {
      // F^-1(p) = - G^-1(p)
      return -erf_inverse_help( sqrt(-2.0*log(p)) );
    }
  else
    {
      // F^-1(p) = G^-1(1-p)
      return erf_inverse_help( sqrt(-2.0*log(1-p)) );
    }
}

double gaussian_distribution(double p, double mu, double sigma)
{
  return sigma*erf_inverse(p)+mu;
}

template<typename scalartype>
bool ASSERT_NON_ZERO(scalartype z) 
{ 
  return bool(std::fabs(z)>1.e-16);
}

template<typename scalartype>
bool ASSERT_NON_ZERO(std::complex<scalartype> z) 
{ 
  return bool(abs(z)>1.e-16);
} 

bool complex_pairs_abs(std::pair<std::complex<double>, int> const& x, 
		       std::pair<std::complex<double>, int> const& y)
{
  return abs(x.first) < abs(y.first);
}

bool complex_less(std::complex<double> const& x, std::complex<double> const& y)
{
return x.real() > y.real();
}

bool complex_less_pairs(std::pair<std::complex<double>, int> const& x, 
			std::pair<std::complex<double>, int> const& y)
{
return x.first.real() < y.first.real();
}



bool pair_less(std::pair<std::complex<double>, std::complex<double> > const& x, std::pair<std::complex<double>, std::complex<double> > const& y)
{
  return abs(x.first) > abs(y.first);
}



bool susceptibility_less(std::complex<double> const& x, std::complex<double> const& y)
{
  return sqrt(square(real(x)-1.) + square(imag(x))) > sqrt(square(real(y)-1.) + square(imag(y)));
}

bool susceptibility_less_pairs(std::pair<std::complex<double>, int> const& x, 
			       std::pair<std::complex<double>, int> const& y)
{
  return sqrt(square(real(x.first)-1.) + square(imag(x.first))) > sqrt(square(real(y.first)-1.) + square(imag(y.first)));
}

bool real_pair_less(std::pair<double, int> const& x,
		    std::pair<double, int> const& y)
{
  return x.first < y.first;
}

/*
double L2_norm(std::vector<double> v1) 
{ 
  double result=0;

  for(size_t i=0; i<v1.size(); i++)
    result += square(v1[i]);

  return result;
}

double L2_norm(std::vector<double> v1, std::vector<double> v2) 
{ 
  assert(v1.size() == v2.size());

  double result=0;

  for(size_t i=0; i<v1.size(); i++)
    result += square(v1[i]-v2[i]);

  return result;
}
*/

/*
bool SAME_VECTOR(std::vector<double> v1, std::vector<double> v2) 
{ 
  return (L2_norm(v1,v2)<1.e-6);
}
*/

template<int DIMENSION>
int column_major_index_2_row_major_index(int /*column_index*/, int* /*grid*/)
{
  throw std::logic_error(__FUNCTION__);
  return -1;
}

template<int DIMENSION>
int row_major_index_2_column_major_index(int /*column_index*/, int* /*grid*/)
{
  throw std::logic_error(__FUNCTION__);
  return -1;
}

template<>
int column_major_index_2_row_major_index<1>(int column_index, int* /*grid*/)
{
  return column_index;
}

template<>
int column_major_index_2_row_major_index<2>(int column_index, int* grid)
{
  int i0 = column_index     % grid[0];
  int i1 = (column_index-i0)/grid[0];

  int row_major_index = grid[1]*i0 + i1;

  return row_major_index;
}

template<>
int column_major_index_2_row_major_index<3>(int column_index, int* grid)
{
  int i0 =  (column_index                       ) % grid[0];
  int i1 = ((column_index-i0           )/grid[0]) % grid[1];
  int i2 = (column_index-i0-i1*grid[0])/(grid[0]*grid[1]);

  int row_major_index = (i0*grid[1] + i1)*grid[2] + i2;

  return row_major_index;
}

template<>
int row_major_index_2_column_major_index<1>(int row_index, int* /*grid*/)
{
  return row_index;
}

template<>
int row_major_index_2_column_major_index<2>(int row_index, int* grid)
{
  int i1 = row_index     %grid[1];
  int i0 = (row_index-i1)/grid[1];

  int column_major_index = i0 + i1*grid[0];

  return column_major_index;
}

template<>
int row_major_index_2_column_major_index<3>(int row_index, int* grid)
{
  int i2 =  (row_index                       ) % grid[2];
  int i1 = ((row_index-i2           )/grid[2]) % grid[1];
  int i0 = (row_index-i2-i1*grid[2])/(grid[1]*grid[2]);

  int column_major_index = i0 + grid[0]*(i1 + grid[1]*i2);

  return column_major_index;
}

bool test_row_major_index_versus_column_major_index()
{
  /*
  cout << __FUNCTION__ << endl;
  
  {
    int* grid = new int[2];
    grid[0] = 7;
    grid[1] = 9;

    for(int l=0; l<grid[0]*grid[1]; l++)
      {
	int row_index = l;
	int col_index = row_major_index_2_column_major_index<2>(row_index, grid);
	assert(row_index == column_major_index_2_row_major_index<2>(col_index, grid));
      }
  }

  {
    int* grid = new int[3];
    grid[0] = 7;
    grid[1] = 9;
    grid[2] = 13;

    for(int l=0; l<grid[0]*grid[1]*grid[2]; l++)
      {
	int row_index = l;
	int col_index = row_major_index_2_column_major_index<3>(row_index, grid);
	assert(row_index == column_major_index_2_row_major_index<3>(col_index, grid));
      }
  }
  */
  return true;
}
