//-*-C++-*-

#ifndef QMC_RESIZEABLE_SQUARE_MATRICES_CPU_H
#define QMC_RESIZEABLE_SQUARE_MATRICES_CPU_H
  
/*! \file   
 *
 *  \author Peter Staar
 *
 *  \brief  Contains a class to represent resizeable square matrices
 */
template<typename scalartype>
class resizeable_square_matrix<scalartype, CPU>
{

public:

  typedef scalartype                                matrix_scalartype;
  typedef resizeable_square_matrix<scalartype, CPU> this_type;

public:

  resizeable_square_matrix();

  resizeable_square_matrix(int current_size);

  resizeable_square_matrix(int current_size,
				    int global_size);

  ~resizeable_square_matrix();

  scalartype& operator()(int i, int j);

  this_type& operator=(this_type& f_other);
  this_type& operator=(this_type f_other);

  int&        get_current_size();
  int         get_global_size();

  void        get_size(std::pair<int,int>& current, std::pair<int,int>& global);

  void resize(int new_global_size);
  void resize_no_copy(int new_current_size);
  void resize_no_copy(int new_current_size, int new_global_size);

  void copy_from(this_type& other_matrix);

  void remove_row_and_column(int i);
  void remove_row_and_column(int index_row, int index_column);

  void swap_row_and_column(int i, int j);
  void cycle_column_forward();
  void cycle_column_backward();
  void copy_row_and_column_into(int i, int j);

  void insert_row_and_column(int index);
  void insert_row_and_column(int index_row, int index_column);

  void insert_row_and_add_column(int index_row);
  void add_row_and_insert_column(int index_column);

  bool difference(this_type& mat);
  void print();

private:

  int                 current_size;
  int                 global_size;

  scalartype*         data;
};

template<typename scalartype>
resizeable_square_matrix<scalartype, CPU>::resizeable_square_matrix():
  current_size(0),
  global_size(64)
{
  data = new scalartype[global_size*global_size];

  memset(data, 0, sizeof(scalartype)*global_size*global_size);
}

template<typename scalartype>
resizeable_square_matrix<scalartype, CPU>::resizeable_square_matrix(int crrnt_size):
  current_size(crrnt_size),
  global_size(64)
{
  assert(current_size>-1);

  if(current_size>global_size)
    global_size = current_size;

  data = new scalartype[global_size*global_size];

  memset(data, 0, sizeof(scalartype)*global_size*global_size);
}

template<typename scalartype>
resizeable_square_matrix<scalartype, CPU>::resizeable_square_matrix(int crrnt_size,
							       int glbl_size):
  current_size(crrnt_size),
  global_size(glbl_size)
{
  assert(current_size>-1);
  assert(global_size >-1);

  if(current_size>global_size)
    global_size = current_size;
 
  data = new scalartype[global_size*global_size];

  memset(data, 0, sizeof(scalartype)*global_size*global_size);
}

template<  typename scalartype>
resizeable_square_matrix<scalartype, CPU>::~resizeable_square_matrix()
{
  delete [] data;
}

template<  typename scalartype>
int& resizeable_square_matrix<scalartype, CPU>::get_current_size()
{
  return current_size;
}

template<  typename scalartype>
int resizeable_square_matrix<scalartype, CPU>::get_global_size()
{
  return global_size;
}

template<  typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::get_size(std::pair<int,int>& current, std::pair<int,int>& global)
{
  current.first  = current_size;
  current.second = current_size;

  global.first  = global_size;
  global.second = global_size;
}

template<  typename scalartype>
inline scalartype& resizeable_square_matrix<scalartype, CPU>::operator()(int i, int j)
{
  assert(i>-1 && j>-1 && i<current_size && j<current_size);
  return data[i + j * global_size];
}

template<typename scalartype>
resizeable_square_matrix<scalartype, CPU>& resizeable_square_matrix<scalartype, CPU>::operator=(resizeable_square_matrix<scalartype, CPU>& M_other)
{
  current_size = M_other.get_current_size();
  global_size  = M_other.get_global_size();
    
  delete [] data;
  data = new scalartype[global_size*global_size];

  return *this;
}

template<typename scalartype>
resizeable_square_matrix<scalartype, CPU>& resizeable_square_matrix<scalartype, CPU>::operator=(resizeable_square_matrix<scalartype, CPU> M_other)
{
  current_size = M_other.get_current_size();
  global_size = M_other.get_global_size();
    
  delete [] data;
  data = new scalartype[global_size*global_size];

  return *this;
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::resize(int new_current_size)
{
  if(new_current_size > global_size)
    {
      int new_global_size = int(1.2*double(new_current_size));

      scalartype* new_data = new scalartype[new_global_size*new_global_size];
      memset(new_data, 0, sizeof(scalartype)*new_global_size*new_global_size);
	
      // copy column by column
      for(int i=0; i<current_size; i++)
	memcpy(&new_data[i*new_global_size], &data[i*global_size], sizeof(scalartype)*current_size);
	
      delete [] data;

      // is this efficient ??? can I not just swap pointers and kill pointer new_data ???
      data = new scalartype[new_global_size*new_global_size];
      memcpy(data, new_data, sizeof(scalartype)*new_global_size*new_global_size);
	
      delete [] new_data;
	
      global_size  = new_global_size;
      current_size = new_current_size;
    }
  else      
    {
      current_size = new_current_size;
    }
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::resize_no_copy(int new_current_size)
{
  if(new_current_size > global_size)
    {
      current_size = new_current_size;
      global_size = int(1.2*double(current_size));

      delete [] data;
      data = new scalartype[global_size*global_size];
    }
  else
    {
      current_size = new_current_size;
    }
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::resize_no_copy(int new_current_size,
							  int new_global_size)
{
  current_size = new_current_size;
  global_size  = new_global_size;

  delete [] data;
  data = new scalartype[global_size*global_size];
    
  memset(data, 0, sizeof(scalartype)*global_size*global_size);
}
 
template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::copy_from(this_type& other_matrix)
{
  current_size = other_matrix.get_current_size();

  if(current_size > 0)
    {
      if(global_size >= other_matrix.get_current_size()) // OK we have space !!
	{
	  if(global_size == other_matrix.get_global_size()) // jack-pot
	    memcpy(data, &(other_matrix(0,0)), sizeof(scalartype)*global_size*global_size);
	  else
	    { // do column by column copy
	      for(int i=0; i<current_size; i++)
		memcpy(&data[i*global_size], &other_matrix(0,i), sizeof(scalartype)*current_size);
	    }
	}
      else // you lose a lot of time here !! avoid at all costs !!
	{
	  global_size = other_matrix.get_global_size();

	  delete [] data;
	  data = new scalartype[global_size*global_size];
	  memcpy(data, &(other_matrix(0,0)), sizeof(scalartype)*global_size*global_size);
	}
    }

  assert(current_size <= global_size);
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::remove_row_and_column(int index)
{
  assert(index>=0);
  assert(index<current_size);

  // copy the remaining columns on the right
  memmove(&data[index*global_size], &data[(index+1)*global_size], sizeof(scalartype)*global_size*(current_size-1-index));

  // copy the remaining rows on the bottom
  for(int i=0; i<current_size-1; i++){
    memmove(&data[index + i*global_size], &data[1 + index + i*global_size], sizeof(scalartype)*(current_size-1-index));
  }

  current_size -= 1;
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::remove_row_and_column(int index_row, int index_column)
{
  assert(index_row>=0 && index_column>=0);
  assert(index_row<current_size && index_column<current_size);

  // copy the remaining columns on the right
  memmove(&data[index_column*global_size], &data[(index_column+1)*global_size], sizeof(scalartype)*global_size*(current_size-1-index_column));

  // copy the remaining rows on the bottom
  for(int i=0; i<current_size-1; i++){
    memmove(&data[index_row + i*global_size], &data[1 + index_row + i*global_size], sizeof(scalartype)*(current_size-1-index_row));
  }

  current_size -= 1;
}


template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::swap_row_and_column(int i, int j)
{
  assert(i>=0 && i<current_size && j>=0 && j<current_size);
    
  //     swap_plan::execute(current_size, &data[i], global_size, &data[j], global_size);
  //     swap_plan::execute(current_size, &data[i*global_size], 1, &data[j*global_size], 1);

    
  scalartype tmp_row, tmp_col;
  for(int l=0; l<current_size; l++)
    {
      tmp_row = data[i + l*global_size];
      data[i + l*global_size] = data[j + l*global_size];
      data[j + l*global_size] = tmp_row;
    }

  for(int l=0; l<current_size; l++)
    {
      tmp_col = data[l + i*global_size];
      data[l + i*global_size] = data[l + j*global_size];
      data[l + j*global_size] = tmp_col;
    }   
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::cycle_column_forward()
{
  scalartype* tmp_column = new scalartype[current_size];

  for(int l=0; l<current_size; l++)
    tmp_column[l] = data[l + (current_size-1)*global_size];

  memmove(&data[(1)*global_size], &data[(0)*global_size], sizeof(scalartype)*global_size*(current_size-1));

  for(int l=0; l<current_size; l++)
    data[l + 0*global_size] = tmp_column[l];

  delete [] tmp_column;
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::cycle_column_backward()
{
  scalartype* tmp_column = new scalartype[current_size];

  for(int l=0; l<current_size; l++)
    tmp_column[l] = data[l + 0*global_size];

  memmove(&data[(0)*global_size], &data[(1)*global_size], sizeof(scalartype)*global_size*(current_size-1));

  for(int l=0; l<current_size; l++)
    data[l + (current_size-1)*global_size] = tmp_column[l];

  delete [] tmp_column;
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::copy_row_and_column_into(int i, int j)
{
  assert(i>=0 && i<current_size && j>=0 && j<current_size);
  // j --> i

  // column
  memcpy(&data[i*global_size], &data[j*global_size], sizeof(scalartype)*(current_size));

  // row 
  for(int l=0; l<current_size; l++)
    data[i + l*global_size] = data[j + l*global_size];
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::insert_row_and_column(int index)
{
  assert(index>=0 && index<current_size);

  resize(current_size+1);

  memmove(&data[(index+1)*global_size], &data[(index)*global_size], sizeof(scalartype)*global_size*(current_size-1-index));

  for(int i=0; i<current_size; i++){
    memmove(&data[index+1 + i*global_size], &data[index + i*global_size], sizeof(scalartype)*(current_size-index));
  }

  for(int l=0; l<current_size; l++){
    data[index + l*global_size] = 0.;
    data[l + index*global_size] = 0.;
  }
}

template<  typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::insert_row_and_column(int index_row, int index_column)
{
  assert(index_row>=0 && index_column>=0);
  assert(index_row<current_size && index_column<current_size);

  resize(current_size+1);

  // copy the remaining columns on the right
  memmove(&data[(index_column+1)*global_size], &data[(index_column)*global_size], sizeof(scalartype)*global_size*(current_size-1-index_column));

  // copy the remaining rows on the bottom
  for(int i=0; i<current_size; i++){
    memmove(&data[(index_row+1) + i*global_size], &data[index_row + i*global_size], sizeof(scalartype)*(current_size-index_row));
  }

  for(int l=0; l<current_size; l++){
    data[index_row + l*global_size] = 0.;
    data[l + index_column*global_size] = 0.;
  }
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::insert_row_and_add_column(int index_row)
{
  assert(index_row>=0 && index_row<current_size);
        
  resize(current_size+1);
        
  // copy the remaining rows on the bottom
  for(int i=0; i<current_size; i++){
    memmove(&data[(index_row+1) + i*global_size], &data[index_row + i*global_size], sizeof(scalartype)*(current_size-index_row));   
  }
        
  for(int l=0; l<current_size; l++){
    data[index_row + l*global_size] = 0.;
    data[l + (current_size-1)*global_size] = 0.;
  }
}
    
template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::add_row_and_insert_column(int index_column)
{
  assert(index_column>=0 && index_column<current_size);
        
  resize(current_size+1);
        
  // copy the remaining columns on the right
  memmove(&data[(index_column+1)*global_size], &data[(index_column)*global_size], sizeof(scalartype)*global_size*(current_size-1-index_column));
        
  for(int l=0; l<current_size; l++){
    data[(current_size-1) + l*global_size] = 0.;
    data[l + index_column*global_size] = 0.;
  }
}
    

template<typename scalartype>
bool resizeable_square_matrix<scalartype, CPU>::difference(this_type& mat)
{
  if(current_size != mat.get_current_size())
    throw std::logic_error(__FUNCTION__);

  cout << scientific;
  cout.precision(8);

  bool OK = true ;
  double max = 0;

  //    int I,J;

  for(int i=0; i<current_size; i++){
    for(int j=0; j<current_size; j++){

      if(max < fabs(mat(i,j)-operator()(i,j)) ){
	// 	  I = i;
	// 	  J = j;
	max = fabs(mat(i,j)-operator()(i,j));
      }

    }
  }

  //if(current_size > 0)
  //max = max/fabs(mat(I,J));

  //cout << "\t\t Max Diff : " << max << endl;

  //     if(max > 1.e-12){
  cout << "\t\t Max Diff : " << max << endl;
  //     }
  //     else
  //       cout << "\t\t\t\tOK\n";

  //     if(max > 1.e-6){
  //       OK = false;
  //       //throw std::logic_error(__FUNCTION__);
  //     }
  /*
    if(!OK)
    {
    for(int i=0; i<current_size; i++){
    for(int j=0; j<current_size; j++){
	    
    if(fabs(mat(i,j)-operator()(i,j)) > 1.e-6)
    cout << "\t" << i<<";" << j;
    else
    cout << "\t --- ";
    }
    cout << endl;
    }
    cout << endl;
    }
  */

  return OK;
}

template<typename scalartype>
void resizeable_square_matrix<scalartype, CPU>::print()
{
    
  std::stringstream ss;
  ss << scientific;
  ss.precision(6);
  ss.width(6);
 
  //     ss << "\n\n";

  /*for(int i=0; i<configuration.size(); i++)
    ss << "\t" << "==============";

    ss << "\n";

    for(int i=0; i<configuration.size(); i++)
    ss << "\t" << configuration[i].get_tau();

    ss << "\n";

    for(int i=0; i<configuration.size(); i++)
    ss << "\t\t" << configuration[i].get_site();

    ss << "\n";

    for(int i=0; i<configuration.size(); i++)
    ss << "\t\t" << configuration[i].get_spin();

    ss << "\n";
  */

  ss << "current-size : " << current_size << "\t global-size : " << global_size << "\n";
  for(int i=0; i<current_size; i++)
    ss << "\t" << "--------------";

  ss << "\n";

  for(int i=0; i<current_size; i++){
    for(int j=0; j<current_size; j++){
      ss << "\t" << this->operator()(i,j);
    }
    ss << endl;
  }
  ss << endl;

  for(int i=0; i<current_size; i++)
    ss << "\t" << "==============";

  ss << endl<<endl;

  cout << ss.str(); 
}

#endif
