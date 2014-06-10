//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file   
 *
 *  Contains a class to represent straightforward matrices
 */

#ifndef QMC_RAW_MATRICES_H
#define QMC_RAW_MATRICES_H

namespace QMC {

  template<typename scalartype>
  class raw_matrix
  {
  public:
    
    raw_matrix(int row_dim, int col_dim);
    ~raw_matrix();
    
    int rows();
    int cols();

    scalartype& operator()(int i, int j);

    void print();

  private:

    int n_rows;
    int n_cols;

    scalartype* data;
  };

  template<typename scalartype>
  raw_matrix<scalartype>::raw_matrix(int row_dim, int col_dim):
    n_rows(row_dim),
    n_cols(col_dim)
  {
    data = new scalartype[row_dim*col_dim];
    memset(data, 0, sizeof(scalartype)*row_dim*col_dim);
  }

  template<typename scalartype>
  raw_matrix<scalartype>::~raw_matrix()
  {
    delete [] data;
  }
  
  template<typename scalartype>
  int raw_matrix<scalartype>::rows()
  { return n_rows; }

  template<typename scalartype>
  int raw_matrix<scalartype>::cols()
  { return n_cols; }

  template<typename scalartype>
  scalartype& raw_matrix<scalartype>::operator()(int i, int j)
  {
    assert(i>-1 && j>-1 && i<n_rows && j<n_cols);
    return data[i + j*n_rows];
  }

  template<typename scalartype>
  void raw_matrix<scalartype>::print()
  {
    cout << "\t";
    for(int j=0; j<n_cols; j++)
      cout << "\t" << j << "\t";
    cout << "\n\n";
    
    for(int i=0; i<n_rows; i++){
      cout << "\t" << i;
      for(int j=0; j<n_cols; j++){
	cout << "\t" << data[i + j*n_rows];
      }
      cout << "\n";
    }
    cout << "\n\n";
  }

}
#endif
