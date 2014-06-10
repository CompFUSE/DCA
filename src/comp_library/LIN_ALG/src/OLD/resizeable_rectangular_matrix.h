//-*-C++-*-

/*! \file   
 *
 * author Peter Staar
 *
 *  Contains a class to represent resizeable rectangular matrices
 */

#ifndef QMC_RESIZEABLE_RECTANGULAR_MATRICES_H
#define QMC_RESIZEABLE_RECTANGULAR_MATRICES_H

namespace QMC {
  
  template<typename scalartype>
  class resizeable_rectangular_matrix
  {
  public:

    typedef scalartype                                matrix_scalartype;
    typedef resizeable_rectangular_matrix<scalartype> this_type;

  public:

    resizeable_rectangular_matrix();

    resizeable_rectangular_matrix(int rows,
				  int cols);

    resizeable_rectangular_matrix(std::pair<int,int> current_size,
				  std::pair<int,int> global_size);

    ~resizeable_rectangular_matrix();

    scalartype& operator()(int i, int j);

    this_type& operator=(this_type& f_other);
    this_type& operator=(this_type f_other);

    std::pair<int,int>& get_current_size();
    std::pair<int,int>  get_global_size();

    void get_size(std::pair<int,int>& current, std::pair<int,int>& global);

    void resize(std::pair<int,int> new_global_size);

    void resize_rows         (std::pair<int,int> new_global_size);
    void resize_cols         (std::pair<int,int> new_global_size);
    void resize_rows_and_cols(std::pair<int,int> new_global_size);


    void resize_no_copy(std::pair<int,int> new_current_size);
    void resize_no_copy(std::pair<int,int> new_current_size, std::pair<int,int> new_global_size);

    void copy_from           (this_type& other_matrix);
    void copy_from_and_resize(this_type& other_matrix);

    void print();

  private:

    std::pair<int,int> current_size;
    std::pair<int,int> global_size;

    scalartype*         data;
  };

  template<typename scalartype>
  resizeable_rectangular_matrix<scalartype>::resizeable_rectangular_matrix()
  {
    current_size.first  = 0;
    current_size.second = 0;
    
    global_size.first  = 0;
    global_size.second = 0;

    data = new scalartype[global_size.first*global_size.second];

    for(int l=0; l<global_size.first*global_size.second; l++)
      data[l] = scalartype(0);
  }

  template<typename scalartype>
  resizeable_rectangular_matrix<scalartype>::resizeable_rectangular_matrix(int rows,
									   int cols)
  {
    current_size.first  = rows;
    current_size.second = cols;
    
    global_size.first  = rows;
    global_size.second = cols;

    data = new scalartype[global_size.first*global_size.second];

    for(int l=0; l<global_size.first*global_size.second; l++)
      data[l] = scalartype(0);
  }

  template<typename scalartype>
  resizeable_rectangular_matrix<scalartype>::resizeable_rectangular_matrix(std::pair<int,int> crrnt_size,
									   std::pair<int,int> glbl_size):
    current_size(crrnt_size),
    global_size(glbl_size)
  {
    if( current_size.first > 0.8*double(global_size.first) )
      global_size.first = int(1.2*double(current_size.first));

    if( current_size.second > 0.8*double(global_size.second) )
      global_size.second = int(1.2*double(current_size.second));


    data = new scalartype[global_size.first*global_size.second];

    for(int l=0; l<global_size.first*global_size.second; l++)
      data[l] = scalartype(0);
  }



  template<  typename scalartype>
  resizeable_rectangular_matrix<scalartype>::~resizeable_rectangular_matrix()
  {
    delete [] data;
  }

  template<  typename scalartype>
  std::pair<int,int>& resizeable_rectangular_matrix<scalartype>::get_current_size()
  {
    return current_size;
  }

  template<  typename scalartype>
  std::pair<int,int> resizeable_rectangular_matrix<scalartype>::get_global_size()
  {
    return global_size;
  }

  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::get_size(std::pair<int,int>& current, std::pair<int,int>& global)
  {
    current.first = current_size.first;
    current.second = current_size.second;
    global.first = global_size.first;
    global.second = global_size.second;
  }

  template<  typename scalartype>
  inline scalartype& resizeable_rectangular_matrix<scalartype>::operator()(int i, int j)
  {
    assert(i >= 0 && j >= 0 && i < current_size.first && j < current_size.second);
    return data[i + j * global_size.first];
  }

  template<typename scalartype>
  resizeable_rectangular_matrix<scalartype>& resizeable_rectangular_matrix<scalartype>::operator=(resizeable_rectangular_matrix<scalartype>& M_other)
  {
    current_size = M_other.get_current_size();
    global_size = M_other.get_global_size();
    
    delete [] data;
    data = new scalartype[global_size.first*global_size.second];

    return *this;
  }

  template<typename scalartype>
  resizeable_rectangular_matrix<scalartype>& resizeable_rectangular_matrix<scalartype>::operator=(resizeable_rectangular_matrix<scalartype> M_other)
  {
    current_size = M_other.get_current_size();
    global_size = M_other.get_global_size();
    
    delete [] data;
    data = new scalartype[global_size.first*global_size.second];

    return *this;
  }

  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::resize(std::pair<int,int> new_current_size)
  {
    if(new_current_size.first > global_size.first || new_current_size.second > global_size.second)
      {
	if(new_current_size.first > global_size.first && new_current_size.second <= global_size.second)
	  resize_rows(new_current_size);

	if(new_current_size.first <= global_size.first && new_current_size.second > global_size.second)
	  resize_cols(new_current_size);
	
	if(new_current_size.first > global_size.first && new_current_size.second > global_size.second)
	  resize_rows_and_cols(new_current_size);
      }
    else      
      {
	current_size = new_current_size;
      }
  }
  
  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::resize_rows(std::pair<int,int> new_current_size)
  {
    // we are ading rows to the matrix

    assert(new_current_size.first > global_size.first && new_current_size.second <= global_size.second);
    
    std::pair<int,int> new_global_size(int(1.2*double(new_current_size.first)), global_size.second);

    scalartype* new_data = new scalartype[new_global_size.first*new_global_size.second];

    for(int j=0; j<current_size.second; j++){

      memcpy(&new_data[j*new_global_size.first], &data[j*global_size.first], sizeof(scalartype)*current_size.first);

      for(int i=current_size.first; i<new_current_size.first; i++)
	new_data[i+j*new_global_size.first] = scalartype(0);
    }

    {   
      scalartype* tmp_ptr = data;
      data     = new_data;
      new_data = tmp_ptr;
    }

    delete [] new_data;

    global_size  = new_global_size;
    current_size = new_current_size;
  }
  
  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::resize_cols(std::pair<int,int> new_current_size)
  {
    // we are ading cols to the matrix

    assert(new_current_size.first <= global_size.first && new_current_size.second > global_size.second);
  
    std::pair<int,int> new_global_size(global_size.first, int(1.2*double(new_current_size.second)));
  
    scalartype* new_data = new scalartype[new_global_size.first*new_global_size.second];

    memcpy(&new_data[0], &data[0], sizeof(scalartype)*global_size.first*current_size.second);

    for(int i=0; i<current_size.first; i++){
      for(int j=current_size.second; j<new_current_size.second; j++){
	new_data[i+j*new_global_size.first] = scalartype(0);
      }
    }

    {   
      scalartype* tmp_ptr = data;
      data     = new_data;
      new_data = tmp_ptr;
    }

    delete [] new_data;
  
    global_size  = new_global_size;
    current_size = new_current_size;
  }

  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::resize_rows_and_cols(std::pair<int,int> new_current_size)
  {
    assert(new_current_size.first > global_size.first && new_current_size.second > global_size.second);

    std::pair<int,int> new_global_size(int(1.2*double(new_current_size.first)), int(1.2*double(new_current_size.second)));
	
    scalartype* new_data = new scalartype[new_global_size.first*new_global_size.second];

    for(int j=0; j<current_size.second; j++)
      memcpy(&new_data[j*new_global_size.first], &data[j*global_size.first], sizeof(scalartype)*current_size.first);

    for(int i=0; i<new_current_size.first; i++){
      for(int j=0; j<new_current_size.second; j++){
	if(i > current_size.first || j > current_size.second)
	  new_data[i+j*new_global_size.first] = scalartype(0);
      }
    }

    {   
      scalartype* tmp_ptr = data;
      data     = new_data;
      new_data = tmp_ptr;
    }

    delete [] new_data;
			   
    global_size  = new_global_size;
    current_size = new_current_size;
  }

  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::resize_no_copy(std::pair<int,int> new_current_size)
  {
    if(new_current_size.first > global_size.first || new_current_size.second > global_size.second)
      {
	current_size = new_current_size;

	if(new_current_size.first > global_size.first && new_current_size.second <= global_size.second)
	  global_size.first = int(1.2*double(new_current_size.first));

	if(new_current_size.first <= global_size.first && new_current_size.second > global_size.second)
	  global_size.second = int(1.2*double(new_current_size.second));
	
	if(new_current_size.first > global_size.first && new_current_size.second > global_size.second){
	  global_size.first  = int(1.2*double(new_current_size.first));
	  global_size.second = int(1.2*double(new_current_size.second));
	}

	delete [] data;
	data = new scalartype[global_size.first*global_size.second];
      }
    else
      {
	current_size = new_current_size;
      }
  }

  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::resize_no_copy(std::pair<int,int> new_current_size,
								 std::pair<int,int> new_global_size)
  {
    current_size = new_current_size;
    global_size  = new_global_size;

    delete [] data;
    data = new scalartype[global_size.first*global_size.second];
  }

 
  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::copy_from(this_type& other_matrix)
  {
    assert(global_size.first == other_matrix.get_global_size().first);
    assert(global_size.second == other_matrix.get_global_size().second);

    current_size = other_matrix.get_current_size();

    if(current_size.first > 0 && current_size.second > 0)
      memcpy(data, &(other_matrix(0,0)), sizeof(scalartype)*global_size.first*global_size.second);
  }

  template<  typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::copy_from_and_resize(this_type& other_matrix)
  {
    current_size = other_matrix.get_current_size();
    global_size  = other_matrix.get_global_size();

    delete [] data;
    data = new scalartype[global_size.first*global_size.second];

    if(current_size.first > 0 && current_size.second > 0)
      memcpy(data, &(other_matrix(0,0)), sizeof(scalartype)*global_size.first*global_size.second);
  }

  template<typename scalartype>
  void resizeable_rectangular_matrix<scalartype>::print()
  {
    cout << "\t";
    for(int j=0; j<current_size.second; j++)
      cout << "\t" << j << "\t";
    cout << "\n\n";
    
    for(int i=0; i<current_size.first; i++){
      cout << "\t" << i;
      for(int j=0; j<current_size.second; j++){
	cout << "\t" << operator()(i,j);
      }
      cout << "\n";
    }
    cout << "\n\n";
  }

}

#endif
