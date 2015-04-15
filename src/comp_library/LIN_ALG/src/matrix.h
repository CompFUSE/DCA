//-*-C++-*-

#ifndef LIN_ALG_MATRICES_H
#define LIN_ALG_MATRICES_H

namespace LIN_ALG {
 
  template<typename scalartype, device_type device_name>
  class matrix
  {
  public:
	
    const static int BLOCK_SIZE = 32;

    typedef matrix<scalartype, device_name> this_type;

    typedef typename MATRIX_SCALARTYPE<scalartype, device_name>::new_scalartype matrix_scalartype;
	
  public:
    
    matrix();
    matrix(std::string name);
	
    matrix(                  int currrent_size);
    matrix(std::string name, int currrent_size);
	
    matrix(                  int currrent_size, int global_size);
    matrix(std::string name, int currrent_size, int global_size);
	
    matrix(                  std::pair<int, int> currrent_size);
    matrix(std::string name, std::pair<int, int> currrent_size);
	
    matrix(                  std::pair<int, int> currrent_size, std::pair<int, int> global_size);
    matrix(std::string name, std::pair<int, int> currrent_size, std::pair<int, int> global_size);

    template<typename other_scalartype, device_type other_device_name>
    matrix(matrix<other_scalartype, other_device_name>& other_matrix);

    ~matrix();
	
    scalartype& operator()(int i, int j);

    std::string& get_name();
    	
    int& get_thread_id();
    int& get_stream_id();

    void set_thread_and_stream_id(int t_id, int s_id);

    matrix_scalartype*& get_ptr();
    matrix_scalartype*  get_ptr(int i, int j);

    bool is_square();

    int get_number_of_rows();
    int get_number_of_cols();
    int get_leading_dimension();

    std::pair<int, int>& get_current_size();
	
    std::pair<int, int>& get_global_size();
      
    void resize(int                 current_size);
    void resize(std::pair<int, int> currrent_size);
      
    void resize_no_copy(int                 current_size);
    void resize_no_copy(std::pair<int, int> currrent_size);
      
    void add_row_and_col();
    void remove_last_row_and_col();
      
    template<device_type other_device_name>
    void swap(matrix<scalartype, other_device_name>& other_matrix);

    template<device_type other_device_name>
    void copy_from(matrix<scalartype, other_device_name>& other_matrix);

    template<device_type other_device_name>
    void copy_from(matrix<scalartype, other_device_name>& other_matrix, copy_concurrency_type copy_t);

    void print();
    void print_fingerprint();
      
    template<device_type other_device_name>
    scalartype difference(matrix<scalartype, other_device_name>& other_matrix);

    template<device_type other_device_name>
    scalartype difference(matrix<scalartype, other_device_name>& other_matrix, int N);
	
  private:
	
    int                 comply_to_block_size(int size);
    std::pair<int, int> comply_to_block_size(std::pair<int, int> size);

    void resize_rows(std::pair<int,int> new_current_size);
    void resize_cols(std::pair<int,int> new_current_size);
	
    void resize_rows_and_cols(std::pair<int,int> new_current_size);
	
  private:

    std::string         name;
	
    int                 thread_id;
    int                 stream_id;

    std::pair<int, int> current_size;
    std::pair<int, int> global_size;
	
    matrix_scalartype*  data;
  };
    
  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix():
    name("unnamed matrix"),
    thread_id(-1),
    stream_id(-1),
    current_size(0 ,0),
    global_size (64,64),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);
    
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::string str):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(0 ,0),
    global_size (64,64),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);
    
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(int c_s):
    name("unnamed matrix"),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s,c_s),
    global_size (c_s,c_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::string str, int c_s):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s,c_s),
    global_size (c_s,c_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(int c_s, int g_s):
    name("unnamed matrix"),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s,c_s),
    global_size (g_s,g_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    assert(c_s>-1 && g_s>-1);
    assert(c_s<=g_s);
	
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::string str, int c_s, int g_s):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s,c_s),
    global_size (g_s,g_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    assert(c_s>-1 && g_s>-1);
    assert(c_s<=g_s);
	
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::pair<int, int> c_s):
    name("unnamed matrix"),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (c_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    assert(c_s.first>-1 && c_s.second>-1);
	
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::string str, std::pair<int, int> c_s):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (c_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    assert(c_s.first>-1 && c_s.second>-1);
	
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }
    
  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::pair<int, int> c_s,
					  std::pair<int, int> g_s):
    name("unnamed matrix"),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (g_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    assert(c_s.first>-1 && c_s.second>-1);
    assert(g_s.first>-1 && g_s.second>-1);
    assert(c_s.first<=g_s.first && c_s.second<=g_s.second);
	
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::matrix(std::string         str, 
					  std::pair<int, int> c_s,
					  std::pair<int, int> g_s):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (g_s),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    assert(c_s.first>-1 && c_s.second>-1);
    assert(g_s.first>-1 && g_s.second>-1);
    assert(c_s.first<=g_s.first && c_s.second<=g_s.second);
	
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size.first*global_size.second);
  }
    
  template<typename scalartype, device_type device_name>
  template<typename other_scalartype, device_type other_device_name>
  matrix<scalartype, device_name>::matrix(matrix<other_scalartype, other_device_name>& other_matrix):	
    name(other_matrix.get_name()),
    thread_id(-1),
    stream_id(-1),
    current_size(other_matrix.get_current_size()),
    global_size (other_matrix.get_global_size()),
    data(NULL)
  {
    global_size = comply_to_block_size(global_size);

    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);

    COPY_FROM<other_device_name, device_name>::execute(other_matrix.get_ptr(), other_matrix.get_current_size(), other_matrix.get_global_size(),
						       data                  , current_size                   , global_size                   );
  }

  template<typename scalartype, device_type device_name>
  matrix<scalartype, device_name>::~matrix()
  {
    MEMORY_MANAGEMENT<device_name>::deallocate(data);
  }

  template<typename scalartype, device_type device_name>
  int& matrix<scalartype, device_name>::get_thread_id()
  {
    return thread_id;
  }

  template<typename scalartype, device_type device_name>
  int& matrix<scalartype, device_name>::get_stream_id()
  {
    return stream_id;
  }

  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::set_thread_and_stream_id(int t_id, int s_id)
  {
    thread_id = t_id;
    stream_id = s_id;
  }
    
  template<typename scalartype, device_type device_name>
  scalartype& matrix<scalartype, device_name>::operator()(int i, int j)
  {
    assert(i>-1 && i<current_size.first);
    assert(j>-1 && j<current_size.second);
    assert(device_name == LIN_ALG::CPU);

    return data[i+j*global_size.first];
  }
   
  template<typename scalartype, device_type device_name>
  std::string& matrix<scalartype, device_name>::get_name()
  {
    return name;
  }
 
  template<typename scalartype, device_type device_name>
  inline typename matrix<scalartype, device_name>::matrix_scalartype*& matrix<scalartype, device_name>::get_ptr()
  {
    return data;
  } 

  template<typename scalartype, device_type device_name>
  inline typename matrix<scalartype, device_name>::matrix_scalartype* matrix<scalartype, device_name>::get_ptr(int i, int j)
  {
    assert(i>-1 && i<current_size.first);
    assert(j>-1 && j<current_size.second);

    return data+(i+j*global_size.first);
  } 

  template<typename scalartype, device_type device_name>
  bool matrix<scalartype, device_name>::is_square()
  {
    return (current_size.first==current_size.second);
  }

  template<typename scalartype, device_type device_name>
  int matrix<scalartype, device_name>::get_number_of_rows()
  {
    return current_size.first;
  }

  template<typename scalartype, device_type device_name>
  int matrix<scalartype, device_name>::get_number_of_cols()
  {
    return current_size.second;
  }
   
  template<typename scalartype, device_type device_name>
  int matrix<scalartype, device_name>::get_leading_dimension()
  {
    return global_size.first;
  }

  template<typename scalartype, device_type device_name>
  std::pair<int,int>& matrix<scalartype, device_name>::get_current_size()
  {
    return current_size;
  }
    
  template<typename scalartype, device_type device_name>
  std::pair<int,int>& matrix<scalartype, device_name>::get_global_size()
  {
    return global_size;
  }

  template<typename scalartype, device_type device_name>
  inline int matrix<scalartype, device_name>::comply_to_block_size(int size)
  {
    assert(size>-1);

    if(size==0)
      return BLOCK_SIZE;
    else{
      int Nb_blocks = (size+BLOCK_SIZE-1)/BLOCK_SIZE;
      return Nb_blocks*BLOCK_SIZE;
    }
  }

  template<typename scalartype, device_type device_name>
  std::pair<int, int> matrix<scalartype, device_name>::comply_to_block_size(std::pair<int, int> size)
  {
    std::pair<int, int> p = size;

    p.first  = comply_to_block_size(p.first );
    p.second = comply_to_block_size(p.second);

    return p;
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize(int current_size)
  {
    std::pair<int, int> c_s(current_size, current_size);
    resize(c_s);
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize(std::pair<int, int> new_current_size)
  {
    if(new_current_size.first > global_size.first || new_current_size.second > global_size.second)
      {
	//CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);

	if(new_current_size.first > global_size.first && new_current_size.second <= global_size.second)
	  resize_rows(new_current_size);
	    
	if(new_current_size.first <= global_size.first && new_current_size.second > global_size.second)
	  resize_cols(new_current_size);
	    
	if(new_current_size.first > global_size.first && new_current_size.second > global_size.second)
	  resize_rows_and_cols(new_current_size);

	//CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);
      }
    else
      {
	current_size = new_current_size;
      }	
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize_rows(std::pair<int,int> new_current_size)
  {
    assert(new_current_size.first > global_size.first && new_current_size.second <= global_size.second);
	
    std::pair<int,int> new_global_size;
    new_global_size.first  = comply_to_block_size(new_current_size.first);
    new_global_size.second = global_size.second;
	
    matrix_scalartype* new_data = NULL;
    MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_global_size);
	
    {
      COPY_FROM<device_name, device_name>::execute(data    , current_size, global_size,
						   new_data, current_size, new_global_size);
    }

    {
      matrix_scalartype* tmp_ptr = data;
      data                       = new_data;
      new_data                   = tmp_ptr;
    }
	
    MEMORY_MANAGEMENT<device_name>::deallocate(new_data);
	
    global_size  = new_global_size;
    current_size = new_current_size;
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize_cols(std::pair<int,int> new_current_size)
  {
    assert(new_current_size.first <= global_size.first && new_current_size.second > global_size.second);
	
    std::pair<int,int> new_global_size;
    new_global_size.first  = global_size.first;
    new_global_size.second = comply_to_block_size(new_current_size.second);

    matrix_scalartype* new_data = NULL;
    MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_global_size);
	
    {
      COPY_FROM<device_name, device_name>::execute(data    , current_size,     global_size,
						   new_data, current_size, new_global_size);
    }

    {
      matrix_scalartype* tmp_ptr = data;
      data                       = new_data;
      new_data                   = tmp_ptr;
    }
    
    MEMORY_MANAGEMENT<device_name>::deallocate(new_data);    
	
    global_size  = new_global_size;
    current_size = new_current_size;
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize_rows_and_cols(std::pair<int,int> new_current_size)
  {
    assert(new_current_size.first > global_size.first && new_current_size.second > global_size.second);
	
    std::pair<int,int> new_global_size;
    new_global_size.first  = comply_to_block_size(new_current_size.first);
    new_global_size.second = comply_to_block_size(new_current_size.second);
	
    matrix_scalartype* new_data = NULL;
    MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_global_size);
	
    {
      COPY_FROM<device_name, device_name>::execute(data    , current_size,     global_size,
						   new_data, current_size, new_global_size);
    }

    {
      matrix_scalartype* tmp_ptr = data;
      data                       = new_data;
      new_data                   = tmp_ptr;
    }
	
    MEMORY_MANAGEMENT<device_name>::deallocate(new_data);
	
    global_size  = new_global_size;
    current_size = new_current_size;
  }

  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize_no_copy(int new_current_size)
  {
    std::pair<int, int> tmp_current_size(new_current_size, new_current_size);
    resize_no_copy(tmp_current_size);
  }
  
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::resize_no_copy(std::pair<int,int> new_current_size)
  {
    if(new_current_size.first > global_size.first || new_current_size.second > global_size.second)
      {
	//CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);
	
	current_size = new_current_size;
	global_size  = comply_to_block_size(new_current_size);

	MEMORY_MANAGEMENT<device_name>::deallocate(data);  
	MEMORY_MANAGEMENT<device_name>::  allocate(data, global_size);

	//CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);
      }
    else
      {
	current_size = new_current_size;
      }
  }

  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::add_row_and_col()
  {
    if(current_size.first<global_size.first && current_size.second<global_size.second)
      {
	current_size.first++;
	current_size.second++;
      }
    else
      {
	std::pair<int, int> new_current_size(current_size.first+1, current_size.second+1);
	this->resize(new_current_size);
      }

    assert(current_size.first >-1 && current_size.first <=global_size.first);
    assert(current_size.second>-1 && current_size.second<=global_size.second);
  }

  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::remove_last_row_and_col()
  {
    current_size.first--;
    current_size.second--;

    assert(current_size.first >-1 && current_size.first <=global_size.first);
    assert(current_size.second>-1 && current_size.second<=global_size.second);
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  void matrix<scalartype, device_name>::swap(matrix<scalartype, other_device_name>& other_matrix)
  {
    if(device_name != other_device_name)
      throw std::logic_error(__FUNCTION__);

    std::pair<int, int> this_c_s  = this->       get_current_size();
    std::pair<int, int> other_c_s = other_matrix.get_current_size();

    std::pair<int, int> this_g_s  = this->       get_global_size();
    std::pair<int, int> other_g_s = other_matrix.get_global_size();

    scalartype* this_ptr  = this->       get_ptr();
    scalartype* other_ptr = other_matrix.get_ptr();
    
    this->get_ptr()          = other_ptr;
    this->get_global_size()  = other_g_s;
    this->get_current_size() = other_c_s;

    other_matrix.get_ptr()          = this_ptr; 
    other_matrix.get_global_size()  = this_g_s;
    other_matrix.get_current_size() = this_c_s;
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  void matrix<scalartype, device_name>::copy_from(matrix<scalartype, other_device_name>& other_matrix)
  {
    this->resize(other_matrix.get_current_size());
    
    COPY_FROM<other_device_name, device_name>::execute(other_matrix.get_ptr(), other_matrix.get_current_size(), other_matrix.get_global_size(),
						       this->get_ptr(),        this->get_current_size(),        this->get_global_size());
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  void matrix<scalartype, device_name>::copy_from(matrix<scalartype, other_device_name>& other_matrix, copy_concurrency_type copy_t)
  {
    const static device_type device_t = LIN_ALG::CUBLAS_DEVICE_NAME<other_device_name, device_name>::device_t;

    assert(thread_id>-1 and stream_id>-1);

    this->resize(other_matrix.get_current_size());

    switch(copy_t)
      {
      case SYNCHRONOUS:
	COPY_FROM<other_device_name, device_name>::execute(other_matrix.get_ptr(), other_matrix.get_current_size(), other_matrix.get_global_size(),
							   this->get_ptr(),        this->get_current_size(),        this->get_global_size());
	break;

      case ASYNCHRONOUS:
	CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);

	COPY_FROM<other_device_name, device_name>::execute(other_matrix.get_ptr(), other_matrix.get_current_size(), other_matrix.get_global_size(),
							   this->get_ptr(),        this->get_current_size(),        this->get_global_size(),
							   thread_id, stream_id);
	
	CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
	break;

      default:
	throw std::logic_error(__FUNCTION__);
      }
  }


  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  scalartype matrix<scalartype, device_name>::difference(matrix<scalartype, other_device_name>& other_matrix)
  {	
    if(this->get_current_size() != other_matrix.get_current_size())
      {
	throw std::logic_error("different matrix size");
      }

    matrix<scalartype, CPU> cp_this (current_size);
    matrix<scalartype, CPU> cp_other(current_size);

    COPY_FROM<device_name,       CPU>::execute(this  ->get_ptr(), this  ->get_current_size(), this  ->get_global_size(),
					       cp_this.get_ptr(), cp_this.get_current_size(), cp_this.get_global_size());
    
    COPY_FROM<other_device_name, CPU>::execute(other_matrix.get_ptr(), other_matrix.get_current_size(), other_matrix.get_global_size(),
					       cp_other    .get_ptr(), cp_other    .get_current_size(), cp_other    .get_global_size()); 

    scalartype max_dif=0;
    for(int i=0; i<current_size.first; ++i){
      for(int j=0; j<current_size.second; ++j){
	/*
	  if( fabs(cp_this(i,j)-cp_other(i,j)) < 1.e-6)
	  std::cout << "\t" << 0.;
	  else
	  std::cout << "\t" << cp_this(i,j)-cp_other(i,j);
	*/

	if( fabs(cp_this(i,j)-cp_other(i,j)) > max_dif)
	  max_dif = fabs(cp_this(i,j)-cp_other(i,j));
      }
      //std::cout << "\n";
    }
    //std::cout << "\n";

    if(fabs(max_dif)<1.e-8)
      {
	//std::cout << "\t\t Max Diff : OK " << endl;
      }
    else
      {
	//std::cout << "\t\t Max Diff : " <<  max_dif << endl;
    
	//throw std::logic_error(__FUNCTION__);
	
	/*
	this->print();
	
	other_matrix.print();
	
     	std::cout << "\n\n";

	for(int i=0; i<current_size.first; ++i){
	    for(int j=0; j<current_size.second; ++j){
		
		if( fabs(cp_this(i,j)-cp_other(i,j)) < 1.e-6)
		    std::cout << "\t" << 0.;
		else
		    std::cout << "\t" << cp_this(i,j)-cp_other(i,j);
	    }
	    std::cout << "\n";
	}
	std::cout << "\n";
	*/
	
	if(fabs(max_dif)>1.e-3)
	    throw std::logic_error(__FUNCTION__);
	
    }

    return max_dif;
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  scalartype matrix<scalartype, device_name>::difference(matrix<scalartype, other_device_name>& other_matrix, int N)
  {	
    if(this->get_current_size() != other_matrix.get_current_size())
      {
	throw std::logic_error("different matrix size");
      }

    matrix<scalartype, CPU> cp_this (current_size);
    matrix<scalartype, CPU> cp_other(current_size);

    COPY_FROM<device_name,       CPU>::execute(this  ->get_ptr(), this  ->get_current_size(), this  ->get_global_size(),
					       cp_this.get_ptr(), cp_this.get_current_size(), cp_this.get_global_size());
    
    COPY_FROM<other_device_name, CPU>::execute(other_matrix.get_ptr(), other_matrix.get_current_size(), other_matrix.get_global_size(),
					       cp_other    .get_ptr(), cp_other    .get_current_size(), cp_other    .get_global_size()); 

    scalartype max_dif=0;
    for(int i=0; i<N; ++i){
      for(int j=0; j<N; ++j){
	if( fabs(cp_this(i,j)-cp_other(i,j)) > max_dif)
	  max_dif = fabs(cp_this(i,j)-cp_other(i,j));
      }
    }

    if(fabs(max_dif)<1.e-8)
      {
	  std::cout << "\t\t Max Diff : OK " << std::endl;
      }
    else
      {
     	std::cout << "\n\n";

	for(int i=0; i<N; ++i){
	    for(int j=0; j<N; ++j){
		
		if( fabs(cp_this(i,j)-cp_other(i,j)) < 1.e-6)
		    std::cout << "\t" << 0.;
		else
		    std::cout << "\t" << cp_this(i,j)-cp_other(i,j);
	    }
	    std::cout << "\n";
	}
	std::cout << "\n";

	/*
	if(fabs(max_dif)>1.e-12)
	    throw std::logic_error(__FUNCTION__);
	*/
    }

    return max_dif;
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::print()
  {
    MEMORY_MANAGEMENT<device_name>::print(data, current_size, global_size);
  }
    
  template<typename scalartype, device_type device_name>
  void matrix<scalartype, device_name>::print_fingerprint()
  {
    std::stringstream ss;

    ss << "\n\n";
    ss << "\t  name : " << name << "\n";
    ss << "\t  current-size : " << current_size.first << ", " << current_size.second << "\n";
    ss << "\t  global-size  : " << global_size.first  << ", " << global_size.second << "\n";
    ss << "\t  memory-size  : " << global_size.first*global_size.second*sizeof(scalartype)*1.e-6 << "(Mbytes)\n"; 
    ss << "\n\n";    

    std::cout << ss.str();
  }
    
}

#endif
