//-*-C++-*-

/*
 * eigensystem_plan.h
 *      Author: bart ydens, Peter Staar
 */


#ifndef BLOCK_EIGENSYSTEM_PLAN_H_
#define BLOCK_EIGENSYSTEM_PLAN_H_

template<typename real_scalartype, matrix_form_type matrix>
class block_eigensystem_plan
{};

template<typename scalartype>
class block_eigensystem_plan< scalartype, GENERAL>
{

public:
  
  block_eigensystem_plan(int* s);
  ~block_eigensystem_plan();
  
  template<typename matrix_type, typename scalartype_2, class domain>
  void execute_plan(function<matrix_type              , dmn_2<domain, domain> >& Matrix,
		    function<matrix_type              , dmn_2<domain, domain> >& V,
		    function<std::vector<scalartype_2>, dmn_2<domain, domain> >& lambda);
  
private:
  int* sizes;
  eigensystem_plan<scalartype, GENERAL> eigen_p;
};

template<typename scalartype>
block_eigensystem_plan< scalartype, GENERAL>::block_eigensystem_plan(int* s):
  sizes(s),
  eigen_p()
{}

template< typename scalartype>
block_eigensystem_plan<scalartype, GENERAL>::~block_eigensystem_plan()
{}

template< typename scalartype>
template<typename matrix_type, typename scalartype_2, class domain>
void block_eigensystem_plan<scalartype, GENERAL>::execute_plan(function<matrix_type              , dmn_2<domain, domain> >& Matrix,
							       function<matrix_type              , dmn_2<domain, domain> >& V,
							       function<std::vector<scalartype_2>, dmn_2<domain, domain> >& lambda)  
{
  std::pair<int,int> Matrix_current_size;
  std::pair<int,int> Matrix_global_size;
  std::pair<int,int> V_current_size;
  std::pair<int,int> V_global_size;
  
  for(int j = 0 ; j< domain::dmn_size();j++)
    {
      for(int i = 0 ; i< domain::dmn_size();i++)
	{
	  Matrix(i,j).get_size(Matrix_current_size, Matrix_global_size);
	  if(Matrix_current_size != std::pair<int,int>(0,0))
	    {
	      V(i,j).get_size(V_current_size, V_global_size);
		
	      assert(i==j); // not implemented
	      assert(Matrix_current_size.first == Matrix_current_size.second);
	      assert(Matrix_current_size.first == V_current_size.first);

	      int n = Matrix_current_size.first;
	      eigen_p.set_N(n);

	      memcpy(eigen_p.A, &Matrix(i,j)(0,0), sizeof(scalartype)*Matrix_global_size.first*Matrix_global_size.second);
		
	      eigen_p.execute_plan();
		
	      memcpy(&lambda(i,j)[0], &eigen_p.W[0], sizeof(scalartype)*Matrix_current_size.first);
		
	      memcpy(&V(i,j)(0,0), &eigen_p.VR[0], sizeof(scalartype)*V_global_size.first*V_global_size.second);
	    }
	}
    }
}

template< typename scalartype>
class block_eigensystem_plan< scalartype, HERMITIAN>
{

public:
    
  block_eigensystem_plan(int* s);
  ~block_eigensystem_plan();
  
  template<typename matrix_type, typename scalartype_2, class domain>
  void execute_plan(function<matrix_type              , dmn_2<domain, domain> >& Matrix,
		    function<matrix_type              , dmn_2<domain, domain> >& V,
		    function<std::vector<scalartype_2>, dmn_2<domain, domain> >& lambda);

private:
  int* sizes;
  eigensystem_plan<scalartype, HERMITIAN> eigen_p;
};
  
template<typename scalartype>
block_eigensystem_plan<scalartype, HERMITIAN>::block_eigensystem_plan(int* s):
  sizes(s),
  eigen_p()
{}

template<typename scalartype>
block_eigensystem_plan<scalartype, HERMITIAN>::~block_eigensystem_plan()
{}


template< typename scalartype>
template<typename matrix_type, typename scalartype_2, class domain>
void block_eigensystem_plan< scalartype, HERMITIAN>::execute_plan(function<matrix_type              , dmn_2<domain, domain> >& Matrix,
								  function<matrix_type              , dmn_2<domain, domain> >& V,
								  function<std::vector<scalartype_2>, dmn_2<domain, domain> >& lambda)
{
  std::pair<int,int> Matrix_current_size;
  std::pair<int,int> Matrix_global_size;
  std::pair<int,int> V_current_size;
  std::pair<int,int> V_global_size;
    
  for(int j = 0 ; j< domain::dmn_size();j++)
    {
      for(int i = 0 ; i< domain::dmn_size();i++)
	{
	  Matrix(i,j).get_size(Matrix_current_size, Matrix_global_size);
	  if(Matrix_current_size != std::pair<int,int>(0,0))
	    {
	      V(i,j).get_size(V_current_size, V_global_size);
		
	      assert(i==j); // not implemented
	      assert(Matrix_current_size.first == Matrix_current_size.second);
	      assert(Matrix_current_size.first == V_current_size.first);
		
	      int n = Matrix_current_size.first;
	      eigen_p.set_N(n);
		
	      memcpy(eigen_p.Matrix, &Matrix(i,j)(0,0), sizeof(scalartype)*Matrix_global_size.first*Matrix_global_size.second);
		
	      eigen_p.execute_plan();
		
	      memcpy(&lambda(i,j)[0], &eigen_p.eigenvalues[0], sizeof(scalartype)*Matrix_current_size.first);
		
	      memcpy(&V(i,j)(0,0), &eigen_p.Matrix[0], sizeof(scalartype)*V_global_size.first*V_global_size.second);
	    }
	}
    }
}

template<typename scalartype>
class block_eigensystem_plan< std::complex<scalartype>, HERMITIAN>
{

public:
    
  block_eigensystem_plan(int* s);
  ~block_eigensystem_plan();
  
  template<typename matrix_type, class domain>
  void execute_plan(function<matrix_type            , dmn_2<domain, domain> >& Matrix,
		    function<matrix_type            , dmn_2<domain, domain> >& V,
		    function<std::vector<scalartype>, dmn_2<domain, domain> >& lambda);

private:
  int* sizes;
  eigensystem_plan<std::complex<scalartype>, HERMITIAN> eigen_p;
};
  
template< typename scalartype>
block_eigensystem_plan<std::complex<scalartype>, HERMITIAN>::block_eigensystem_plan(int* s):
  sizes(s),
  eigen_p(1)
{}

template<typename scalartype>
block_eigensystem_plan<std::complex<scalartype>, HERMITIAN>::~block_eigensystem_plan()
{}


template<typename scalartype>
template<typename matrix_type, class domain>
void block_eigensystem_plan<std::complex<scalartype>, HERMITIAN>::execute_plan(function<matrix_type            , dmn_2<domain, domain> >& Matrix,
									       function<matrix_type            , dmn_2<domain, domain> >& V,
									       function<std::vector<scalartype>, dmn_2<domain, domain> >& lambda)
{
  std::pair<int,int> Matrix_current_size;
  std::pair<int,int> Matrix_global_size;
  std::pair<int,int> V_current_size;
  std::pair<int,int> V_global_size;
    
  for(int j = 0 ; j< domain::dmn_size();j++)
    {
      for(int i = 0 ; i< domain::dmn_size();i++)
	{
	  Matrix(i,j).get_size(Matrix_current_size, Matrix_global_size);
	  if(Matrix_current_size != std::pair<int,int>(0,0))
	    {
	      V(i,j).get_size(V_current_size, V_global_size);
		
	      assert(i==j); // not implemented
	      assert(Matrix_current_size.first == Matrix_current_size.second);
	      assert(Matrix_current_size.first == V_current_size.first);
		
	      int n = Matrix_current_size.first;
	      eigen_p.set_N(n);
		
	      memcpy(eigen_p.Matrix, &Matrix(i,j)(0,0), sizeof(std::complex<scalartype>)*Matrix_global_size.first*Matrix_global_size.second);
		
	      eigen_p.execute_plan();
		
	      memcpy(&lambda(i,j)[0], &eigen_p.eigenvalues[0], sizeof(scalartype)*Matrix_current_size.first);
		
	      memcpy(&V(i,j)(0,0), &eigen_p.Matrix[0], sizeof(std::complex<scalartype>)*V_global_size.first*V_global_size.second); 
	    }
	}
    }
}
#endif
