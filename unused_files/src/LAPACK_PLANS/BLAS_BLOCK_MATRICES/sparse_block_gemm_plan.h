//-*-C++-*-

/*
 * sparse_block_matrix_sparse_block_matrix_multiplication_plan.h
 *
 * only applicable to diagonal sparse matices. NO 2 BLOCKS ON SAME ROW OR COLUMN !!!
 *
 *      Author: bart ydens, peter staar
 */

#ifndef SPARSE_BLOCK_MATRIX_MATRIX_PLAN_H_
#define SPARSE_BLOCK_MATRIX_MATRIX_PLAN_H_

template<typename scalartype, class domain>
class sparse_block_gemm_plan
{
  
public:
  
  sparse_block_gemm_plan();
  ~sparse_block_gemm_plan();
  
  template<class matrix_type>
  void execute_plan(sparse_block_matrix<matrix_type, domain>& Matrix_left,
		    sparse_block_matrix<matrix_type, domain>& Matrix_right,
		    sparse_block_matrix<matrix_type, domain>& Matrix_result);

  template<class matrix_type, class vector_type>
  void execute_plan(sparse_block_matrix<matrix_type, domain>& Matrix_left,
		    diagonal_matrix    <vector_type, domain>& Vector_right,
		    sparse_block_matrix<matrix_type, domain>& Matrix_result);

  template<class matrix_type, class vector_type>
  void execute_plan(diagonal_matrix    <vector_type, domain>& Vector_left,
		    sparse_block_matrix<matrix_type, domain>& Matrix_right,
		    sparse_block_matrix<matrix_type, domain>& Matrix_result);

  template<class matrix_type, class vector_type>
  void execute_plan(sparse_block_matrix<matrix_type, domain>& Matrix_left,
		     diagonal_matrix    <vector_type, domain>& Vector_right,
		     sparse_block_matrix<matrix_type, domain>& Matrix_result,
		     std::vector<int> non_zero);
  
private:
  void assert_sizes(std::pair<int,int> m1,
		    std::pair<int,int> m2,
		    std::pair<int,int> m3);
    
private:
  gemm_plan<scalartype> gemm_p;
};

template<typename scalartype, class domain>
sparse_block_gemm_plan<scalartype,domain>::sparse_block_gemm_plan():
  gemm_p()
{}

template<typename scalartype, class domain>
sparse_block_gemm_plan<scalartype,domain>::~sparse_block_gemm_plan()
{}
  
template<typename scalartype, class domain>
template<typename matrix_type>
void sparse_block_gemm_plan<scalartype, domain>::execute_plan(sparse_block_matrix<matrix_type, domain>& Matrix_left,
							      sparse_block_matrix<matrix_type, domain>& Matrix_right,
							      sparse_block_matrix<matrix_type, domain>& Matrix_result)
{
  gemm_p.alpha = 1;
  gemm_p.beta  = 0;

  gemm_p.TRANSA = 'N';
  gemm_p.TRANSB = 'N';

  Matrix_result.blocks.clear();
  int k = 0;
  for(int i = 0; i< (int)   Matrix_left.blocks.size(); i++){
    for(int j = 0; j< (int)   Matrix_right.blocks.size(); j++){
      if(Matrix_left.blocks[i].second == Matrix_right.blocks[j].first){
	Matrix_result.add_new_empty(std::pair<int,int>(Matrix_left.blocks[i].first,Matrix_right.blocks[j].second), 
				    std::pair<int,int>(Matrix_left.values[i].get_current_size().first,Matrix_right.values[j].get_current_size().second) );

	gemm_p.M = Matrix_left. values[i].get_current_size().first;
	gemm_p.N = Matrix_right.values[j].get_current_size().second;
	gemm_p.K = Matrix_left. values[i].get_current_size().second;

	gemm_p.LDA = Matrix_left  .values[i].get_global_size().first;
	gemm_p.LDB = Matrix_right .values[j].get_global_size().first;
	gemm_p.LDC = Matrix_result.values[k].get_global_size().first;
	
	gemm_p.A = &Matrix_left  .values[i](0,0);
	gemm_p.B = &Matrix_right .values[j](0,0);
	gemm_p.C = &Matrix_result.values[k](0,0);

	gemm_p.execute_plan();
	k++;
      }
    }
  }
}

template<typename scalartype, class domain>
template<typename matrix_type, typename vector_type>
void sparse_block_gemm_plan<scalartype, domain>::execute_plan(sparse_block_matrix<matrix_type, domain>& Matrix_left,
							      diagonal_matrix    <vector_type, domain>& Vector_right,
							      sparse_block_matrix<matrix_type, domain>& Matrix_result)

{
  Matrix_result.copy(Matrix_left);

  for(int i = 0; i< (int) Matrix_left.blocks.size(); i++){
    int j = Matrix_left.blocks[i].second;
    for(int b = 0 ;b< Matrix_left(i).get_current_size().second;b++)
	BLAS::SCAL(Matrix_left(i).get_current_size().first, 
		   Vector_right(j)[b], 
		   &Matrix_result(i)(0,b), 
		   1);
  }

}

template<typename scalartype, class domain>
template<typename matrix_type, typename vector_type>
void sparse_block_gemm_plan<scalartype, domain>::execute_plan(diagonal_matrix    <vector_type, domain>& Vector_left,
							      sparse_block_matrix<matrix_type, domain>& Matrix_right,
							      sparse_block_matrix<matrix_type, domain>& Matrix_result)
{
  Matrix_result.copy(Matrix_right);

  for(int i = 0; i< (int) Matrix_right.blocks.size(); i++){
    int j = Matrix_right.blocks[i].first;
    for(int a = 0 ; a< Matrix_right(i).get_current_size().first;a++)
	BLAS::SCAL(Matrix_right(i).get_current_size().second, 
		   Vector_left(j)[a], 
		   &Matrix_result(i)(a,0), 
		   Matrix_right(i).get_current_size().first);
  }
}

#endif
