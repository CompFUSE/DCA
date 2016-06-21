//-*-C++-*-

#ifndef BLOCK_MATRIX_MATRIX_PLAN_H_
#define BLOCK_MATRIX_MATRIX_PLAN_H_


/*!
 * block_matrix_block_matrix_multiplication_plan.h
 *
 *      \author bart ydens, peter staar
 */
template<typename scalartype>
class block_gemm_plan
{
  
public:
  
  block_gemm_plan();
  ~block_gemm_plan();
  
  template<class matrix_type, class domain>
  void execute_plan(function<matrix_type, dmn_2<domain, domain> >& Matrix_1,
		    function<matrix_type, dmn_2<domain, domain> >& Matrix_2,
		    function<matrix_type, dmn_2<domain, domain> >& Matrix_result,
		    scalartype a=1, scalartype b=0,
		    char A='N'    , char B='N');
  
  template<typename matrix_type, class domain, typename other_scalartype>
  void execute_plan(function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
		    function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
		    function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2);
  
  template<typename matrix_type, class domain, typename other_scalartype>
  void execute_plan(function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
		    function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
		    function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2);

  template<typename matrix_type, class domain, typename other_scalartype_0, typename other_scalartype_1>
  void execute_plan_resize(function<matrix_type ,dmn_2<domain, domain> >& Matrix_1,
			   function<matrix_type ,dmn_2<domain, domain> >& Matrix_2,
			   function<matrix_type ,dmn_2<domain, domain> >& Matrix_result,
			   other_scalartype_0 a = 1., other_scalartype_1 b = 0.,
			   char A = 'N', char B = 'N');

  template<typename matrix_type, class domain, typename other_scalartype>
  void execute_plan_resize(function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
			   function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
			   function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2);

  template<typename matrix_type, class domain, typename other_scalartype>
  void execute_plan_resize(function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
			   function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
			   function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2);

private:
  void assert_sizes(std::pair<int,int> m1,
		    std::pair<int,int> m2,
		    std::pair<int,int> m3);
    
private:
  gemm_plan<scalartype> gemm_p;
};

template<typename scalartype>
block_gemm_plan<scalartype>::block_gemm_plan():
  gemm_p()
{}

template<typename scalartype>
block_gemm_plan<scalartype>::~block_gemm_plan()
{}
  
template<typename scalartype>
template<typename matrix_type, class domain>
void block_gemm_plan<scalartype>::execute_plan(function<matrix_type ,dmn_2<domain, domain> >& Matrix_1,
					       function<matrix_type ,dmn_2<domain, domain> >& Matrix_2,
					       function<matrix_type ,dmn_2<domain, domain> >& Matrix_result,
					       scalartype a, scalartype b,
					       char A, char B)
{
  std::pair<int,int> Matrix_1_current_size;
  std::pair<int,int> Matrix_2_current_size;
  std::pair<int,int> Matrix_result_current_size;

  std::pair<int,int> Matrix_1_global_size;
  std::pair<int,int> Matrix_2_global_size;
  std::pair<int,int> Matrix_result_global_size;

  gemm_p.alpha = a;
  gemm_p.beta  = b;
  
  gemm_p.TRANSA = A;
  gemm_p.TRANSB = B;
    
  for(int i = 0 ; i< domain::dmn_size();i++)
    {
      for(int j = 0 ; j< domain::dmn_size();j++)
	{
	  for(int k = 0 ; k< domain::dmn_size();k++)
	    {
	      Matrix_1((A == 'N' ? i : k),(A == 'N' ? k : i)).get_size(Matrix_1_current_size, Matrix_1_global_size);
	      Matrix_2((B == 'N' ? k : j),(B == 'N' ? j : k)).get_size(Matrix_2_current_size, Matrix_2_global_size);

	      if(Matrix_1_current_size != std::pair<int,int>(0,0) && Matrix_2_current_size != std::pair<int,int>(0,0))
		{
		  Matrix_result(i,j).get_size(Matrix_result_current_size, Matrix_result_global_size);
		  assert_sizes(Matrix_1_current_size,Matrix_2_current_size,Matrix_result_current_size);
		    
		  gemm_p.M = Matrix_1_current_size.first;
		  gemm_p.N = Matrix_2_current_size.second;
		  gemm_p.K = Matrix_1_current_size.second;

		  gemm_p.LDA = Matrix_1_global_size.first;
		  gemm_p.LDB = Matrix_2_global_size.first;
		  gemm_p.LDC = Matrix_result_global_size.first;

		  gemm_p.A = &Matrix_1((A == 'N' ? i : k),(A == 'N' ? k : i))(0,0);
		  gemm_p.B = &Matrix_2((B == 'N' ? k : j),(B == 'N' ? j : k))(0,0);
		  gemm_p.C = &Matrix_result(i,j)(0,0);
		    
		  gemm_p.execute_plan();
		}
	    }
	}
    }
}
  
template<typename scalartype>
template<typename matrix_type, class domain, typename other_scalartype>
void  block_gemm_plan<scalartype>::execute_plan(function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
						function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
						function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2)
{
  std::pair<int,int> Matrix_1_current_size;
  
  std::pair<int,int> Matrix_1_global_size;

  for(int i = 0 ; i< domain::dmn_size();i++)
    {
      for(int j = 0 ; j< domain::dmn_size();j++)
	{
	  Matrix_1(i,j).get_size(Matrix_1_current_size, Matrix_1_global_size);
	  if(Matrix_1_current_size != std::pair<int,int>(0,0))
	    {
	      for(int b = 0 ;b< Matrix_1_current_size.second;b++)
		{
		  for(int a = 0 ; a< Matrix_1_current_size.first;a++)
		    {
		      Matrix_2(i,j)(a,b) = Matrix_1(i,j)(a,b) * vector(j,j)[b];
		    }
		}
	    }
	}
    }
 }

template<typename scalartype>
template<typename matrix_type, class domain, typename other_scalartype>
void  block_gemm_plan<scalartype>::execute_plan(function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
						function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
						function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2)
{
  std::pair<int,int> Matrix_1_current_size;
  
  std::pair<int,int> Matrix_1_global_size;
  
  for(int i = 0 ; i< domain::dmn_size();i++)
    {
      for(int j = 0 ; j< domain::dmn_size();j++)
	{
	  Matrix_1(i,j).get_size(Matrix_1_current_size, Matrix_1_global_size);
	  if(Matrix_1_current_size != std::pair<int,int>(0,0))
	    {
	      for(int b = 0 ;b< Matrix_1_current_size.second;b++)
		{
		  for(int a = 0 ; a< Matrix_1_current_size.first;a++)
		    {
		      Matrix_2(i,j)(a,b) = vector(i,i)[a] * Matrix_1(i,j)(a,b);
		    }
		}
	    }
	}
    }
}
 
template<typename scalartype>
template<typename matrix_type, class domain, typename other_scalartype_0, typename other_scalartype_1>
void block_gemm_plan<scalartype>::execute_plan_resize(function<matrix_type ,dmn_2<domain, domain> >& Matrix_1,
						      function<matrix_type ,dmn_2<domain, domain> >& Matrix_2,
						      function<matrix_type ,dmn_2<domain, domain> >& Matrix_result,
						      other_scalartype_0 a, other_scalartype_1 b,
						      char A, char B)
{
  static std::pair<int,int> Matrix_1_current_size;
  static std::pair<int,int> Matrix_2_current_size;
  static std::pair<int,int> Matrix_result_current_size;

  static std::pair<int,int> Matrix_1_global_size;
  static std::pair<int,int> Matrix_2_global_size;
  static std::pair<int,int> Matrix_result_global_size;

  if(b != other_scalartype_1(0) && b != other_scalartype_1(1))
    throw std::logic_error(__FUNCTION__);

  gemm_p.alpha = a;
  gemm_p.beta = 1;
  
  gemm_p.TRANSA = A;
  gemm_p.TRANSB = B;
    
  for(int i = 0 ; i< domain::dmn_size();i++)
    {
      for(int j = 0 ; j< domain::dmn_size();j++)
	{
	  if(b == other_scalartype_1(0)){
	    Matrix_result(i,j).get_size(Matrix_result_current_size, Matrix_result_global_size);
	    if(Matrix_result_current_size != std::pair<int,int>(0,0))
	      Matrix_result(i,j).resize_no_copy(std::pair<int,int>(0,0),std::pair<int,int>(0,0));
	  }
	  for(int k = 0 ; k< domain::dmn_size();k++)
	    {
	      Matrix_1((A == 'N' ? i : k),(A == 'N' ? k : i)).get_size(Matrix_1_current_size, Matrix_1_global_size);
	      Matrix_2((B == 'N' ? k : j),(B == 'N' ? j : k)).get_size(Matrix_2_current_size, Matrix_2_global_size);

	      if(Matrix_1_current_size != std::pair<int,int>(0,0) && Matrix_2_current_size != std::pair<int,int>(0,0))
		{
		  Matrix_result(i,j).get_size(Matrix_result_current_size, Matrix_result_global_size);
		  if(Matrix_result_current_size == std::pair<int,int>(0,0))
		    {		   
		      std::pair<int,int> Matrix_result_size;
		      Matrix_result_size.first  = (A == 'N' ? Matrix_1_current_size.first  : Matrix_1_current_size.second);
		      Matrix_result_size.second = (B == 'N' ? Matrix_2_current_size.second : Matrix_2_current_size.first );
	
		      Matrix_result(i,j).resize_no_copy(Matrix_result_size, Matrix_result_size);
		    }
		  Matrix_result(i,j).get_size(Matrix_result_current_size, Matrix_result_global_size);
		  assert_sizes(Matrix_1_current_size,Matrix_2_current_size,Matrix_result_current_size);
		    
		  gemm_p.M = Matrix_1_current_size.first;
		  gemm_p.N = Matrix_2_current_size.second;
		  gemm_p.K = Matrix_1_current_size.second;

		  gemm_p.LDA = Matrix_1_global_size.first;
		  gemm_p.LDB = Matrix_2_global_size.first;
		  gemm_p.LDC = Matrix_result_global_size.first;

		  gemm_p.A = &Matrix_1((A == 'N' ? i : k),(A == 'N' ? k : i))(0,0);
		  gemm_p.B = &Matrix_2((B == 'N' ? k : j),(B == 'N' ? j : k))(0,0);
		  gemm_p.C = &Matrix_result(i,j)(0,0);
		    
		  gemm_p.execute_plan();
		}
	    }
	}
    }
}

template<typename scalartype>
template<typename matrix_type, class domain, typename other_scalartype>
void  block_gemm_plan<scalartype>::execute_plan_resize(function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
						       function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
						       function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2)
{
  static std::pair<int,int> Matrix_1_current_size;
  static std::pair<int,int> Matrix_2_current_size;

  static std::pair<int,int> Matrix_1_global_size;
  static std::pair<int,int> Matrix_2_global_size;

  for(int i = 0 ; i< domain::dmn_size();i++)
    {
      for(int j = 0 ; j< domain::dmn_size();j++)
	{
	  Matrix_1(i,j).get_size(Matrix_1_current_size, Matrix_1_global_size);
	  if(Matrix_1_current_size != std::pair<int,int>(0,0))
	    {
	      Matrix_2(i,j).resize_no_copy(Matrix_1_current_size,Matrix_1_current_size);
	      for(int b = 0 ;b< Matrix_1_current_size.second;b++)
		{
		  for(int a = 0 ; a< Matrix_1_current_size.first;a++)
		    {
		      Matrix_2(i,j)(a,b) = Matrix_1(i,j)(a,b) * vector(j,j)[b];
		    }
		}
	    }
	  else{
	    Matrix_2(i,j).get_size(Matrix_2_current_size, Matrix_2_global_size);
	    if(Matrix_2_current_size != std::pair<int,int>(0,0))
	      Matrix_2(i,j).resize_no_copy(std::pair<int,int>(0,0),std::pair<int,int>(0,0));
	  }
	}
    }
}

template<typename scalartype>
template<typename matrix_type, class domain, typename other_scalartype>
 void  block_gemm_plan<scalartype>::execute_plan_resize(function<std::vector<other_scalartype>, dmn_2<domain, domain> >& vector,
							function<matrix_type                  , dmn_2<domain, domain> >& Matrix_1,
							function<matrix_type                  , dmn_2<domain, domain> >& Matrix_2)
  {
    static std::pair<int,int> Matrix_1_current_size;
    static std::pair<int,int> Matrix_2_current_size;

    static std::pair<int,int> Matrix_1_global_size;
    static std::pair<int,int> Matrix_2_global_size;

    for(int i = 0 ; i< domain::dmn_size();i++)
      {
	for(int j = 0 ; j< domain::dmn_size();j++)
	  {
	    Matrix_1(i,j).get_size(Matrix_1_current_size, Matrix_1_global_size);
	    if(Matrix_1_current_size != std::pair<int,int>(0,0))
	      {
		Matrix_2(i,j).resize_no_copy(Matrix_1_current_size,Matrix_1_current_size);
		for(int b = 0 ;b< Matrix_1_current_size.second;b++)
		  {
		    for(int a = 0 ; a< Matrix_1_current_size.first;a++)
		      {
			Matrix_2(i,j)(a,b) = vector(i,i)[a] * Matrix_1(i,j)(a,b);
		      }
		  }
	      }
	    else{
	      Matrix_2(i,j).get_size(Matrix_2_current_size, Matrix_2_global_size);
	      if(Matrix_2_current_size != std::pair<int,int>(0,0))
		Matrix_2(i,j).resize_no_copy(std::pair<int,int>(0,0),std::pair<int,int>(0,0));
	    }
	  }
      }
  }




						
template<typename scalartype>
void block_gemm_plan<scalartype>::assert_sizes(std::pair<int,int> m1,
					       std::pair<int,int> m2,
					       std::pair<int,int> m3)
{
  assert(m1.first  == m3.first);
  assert(m1.second == m2.first);
  assert(m2.second == m3.second);
}
#endif
