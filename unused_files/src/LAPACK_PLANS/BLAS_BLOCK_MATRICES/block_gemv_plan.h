//-*-C++-*-

/*
 * block_matrix_block_matrix_multiplication_plan.h
 *
 *      Author: bart ydens, Peter Staar
 */

#ifndef BLOCK_MATRIX_VECTOR_PLAN_H_
#define BLOCK_MATRIX_VECTOR_PLAN_H_

template<typename scalartype>
class block_gemv_plan
{
  
public:
    
  block_gemv_plan(int* s);
  ~block_gemv_plan();
    
  template<typename matrix_type, class domain>
  void execute_plan(FUNC_LIB::function<matrix_type ,dmn_2<domain, domain> >& Matrix,
		    scalartype* vector_source,
		    scalartype* vector_target);
    
private:
  int* sizes;
  gemv_plan<scalartype> gemv_p;
};

template<typename scalartype>
block_gemv_plan<scalartype>::block_gemv_plan(int* s):
  sizes(s),
  gemv_p()
{}

template<typename scalartype>
block_gemv_plan<scalartype>::~block_gemv_plan()
{}
  
template<typename scalartype>
template<typename matrix_type, class domain>
void block_gemv_plan<scalartype>::execute_plan(FUNC_LIB::function<matrix_type, dmn_2<domain, domain> >& Matrix,
					       scalartype*                                    vector_source,
					       scalartype*                                    vector_target)
{
  for(int j = 0 ; j< domain::dmn_size();j++){
    for(int i = 0 ; i< domain::dmn_size();i++){
      if(Matrix(i,j).get_current_size() != std::pair<int,int>(0,0)){
	gemv_p.M = Matrix(i,j).get_current_size().first;
	gemv_p.N = Matrix(i,j).get_current_size().second;
	gemv_p.LDA = Matrix(i,j).get_global_size().first;
	gemv_p.matrix = &Matrix(i,j)(0,0);
	int tmp = 0;
	for(int k = 0; k<i;k++)
	  tmp += sizes[k];
	gemv_p.vector_source = &vector_source[tmp];
	gemv_p.vector_target = &vector_target[tmp];
	  
	gemv_p.execute_plan();
      }
    }
  }
}


#endif
