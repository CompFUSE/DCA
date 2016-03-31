//-*-C++-*-

#ifndef LIN_ALG_VECTORS_H
#define LIN_ALG_VECTORS_H

namespace LIN_ALG {

  template<typename scalartype, device_type device_name>
  class vector
  {
  public:

    typedef vector<scalartype, device_name> this_type;

    typedef typename MATRIX_SCALARTYPE<scalartype, device_name>::new_scalartype vector_scalartype;

  public:

    vector();
    vector(std::string name);

    vector(                  int currrent_size);
    vector(std::string name, int currrent_size);

    vector(                  int currrent_size, int global_size);
    vector(std::string name, int currrent_size, int global_size);

    template<typename other_scalartype, device_type other_device_name>
    vector(vector<other_scalartype, other_device_name>& other_vector);

    ~vector();

    scalartype& operator[](int i);

    std::string& get_name();

    int& get_thread_id();
    int& get_stream_id();

    void set_thread_and_stream_id(int t_id, int s_id);

    void erase    (int i);
    void push_back(scalartype val);

    void set(std::vector<scalartype>& other_vector);
    void set(std::vector<scalartype>& other_vector, copy_concurrency_type copy_t);

    template<device_type other_device_name>
    void set(vector<scalartype, other_device_name>& other_vector);

    template<device_type other_device_name>
    void set(vector<scalartype, other_device_name>& other_vector, copy_concurrency_type copy_t);

    vector_scalartype* get_ptr();
    vector_scalartype* get_ptr(int i);

    int& size();
    int& get_current_size();
    int& get_global_size();

    void reserve(int new_current_size);
    void resize (int new_current_size);

    template<device_type other_device_name>
    scalartype difference(vector<scalartype, other_device_name>& other_vector);

    scalartype difference(std::vector<scalartype>& other_vector);

    void print();

  private:

    std::string        name;

    int                thread_id;
    int                stream_id;

    int                current_size;
    int                global_size;

    vector_scalartype* data;
  };

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::vector():
    name("unnamed vector"),
    thread_id(-1),
    stream_id(-1),
    current_size(0),
    global_size (64),
    data(NULL)
  {
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size);
  }

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::vector(std::string str):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(0),
    global_size (64),
    data(NULL)
  {
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size);
  }

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::vector(int c_s):
    name("unnamed vector"),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (c_s),
    data(NULL)
  {
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size);
  }

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::vector(std::string str, int c_s):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (c_s),
    data(NULL)
  {
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size);
  }

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::vector(int c_s, int g_s):
    name("unnamed vector"),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (g_s),
    data(NULL)
  {
    assert(c_s>-1 && g_s>-1);
    assert(c_s<=g_s);

    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size);
  }

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::vector(std::string str, int c_s, int g_s):
    name(str),
    thread_id(-1),
    stream_id(-1),
    current_size(c_s),
    global_size (g_s),
    data(NULL)
  {
    assert(c_s>-1 && g_s>-1);
    assert(c_s<=g_s);

    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    MEMORY_MANAGEMENT<device_name>::set_to_zero(data, global_size);
  }

  template<typename scalartype, device_type device_name>
  template<typename other_scalartype, device_type other_device_name>
  vector<scalartype, device_name>::vector(vector<other_scalartype, other_device_name>& other_vector):
    name(other_vector.get_name()),

    thread_id(-1),
    stream_id(-1),

    current_size(other_vector.get_current_size()),
    global_size (other_vector.get_global_size()),

    data(NULL)
  {
    MEMORY_MANAGEMENT<device_name>::allocate(data, global_size);
    COPY_FROM<other_device_name, device_name>::execute(other_vector.get_ptr(), this->get_ptr(), current_size);
  }

  template<typename scalartype, device_type device_name>
  vector<scalartype, device_name>::~vector()
  {
    MEMORY_MANAGEMENT<device_name>::deallocate(data);
  }

  template<typename scalartype, device_type device_name>
  scalartype& vector<scalartype, device_name>::operator[](int i)
  {
    assert(i>-1 && i<current_size);
    assert(device_name == LIN_ALG::CPU);

    return data[i];
  }

  template<typename scalartype, device_type device_name>
  std::string& vector<scalartype, device_name>::get_name()
  {
    return name;
  }

  template<typename scalartype, device_type device_name>
  int& vector<scalartype, device_name>::get_thread_id()
  {
    return thread_id;
  }

  template<typename scalartype, device_type device_name>
  int& vector<scalartype, device_name>::get_stream_id()
  {
    return stream_id;
  }

  template<typename scalartype, device_type device_name>
  void vector<scalartype, device_name>::set_thread_and_stream_id(int t_id, int s_id)
  {
    thread_id = t_id;
    stream_id = s_id;
  }

  template<typename scalartype, device_type device_name>
  void vector<scalartype, device_name>::set(std::vector<scalartype>& other_vector)
  {
    reserve(other_vector.size());
    COPY_FROM<LIN_ALG::CPU, device_name>::execute(&other_vector[0], this->get_ptr(), current_size);
  }

  template<typename scalartype, device_type device_name>
  void vector<scalartype, device_name>::set(std::vector<scalartype>& other_vector, copy_concurrency_type copy_t)
  {
    const static device_type device_t = LIN_ALG::CUBLAS_DEVICE_NAME<LIN_ALG::CPU, device_name>::device_t;

    assert(thread_id>-1 and stream_id>-1);

    reserve(other_vector.size());

    switch(copy_t)
    {
    case SYNCHRONOUS:
      COPY_FROM<LIN_ALG::CPU, device_name>::execute(&other_vector[0], this->get_ptr(), current_size);
      break;

    case ASYNCHRONOUS:
      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);

      COPY_FROM<LIN_ALG::CPU, device_name>::execute(&other_vector[0], this->get_ptr(), current_size, thread_id, stream_id);

      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  void vector<scalartype, device_name>::set(vector<scalartype, other_device_name>& other_vector)
  {
    reserve(other_vector.size());
    COPY_FROM<other_device_name, device_name>::execute(other_vector.get_ptr(), this->get_ptr(), current_size);
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  void vector<scalartype, device_name>::set(vector<scalartype, other_device_name>& other_vector, copy_concurrency_type copy_t)
  {
    const static device_type device_t = LIN_ALG::CUBLAS_DEVICE_NAME<other_device_name, device_name>::device_t;

    assert(thread_id>-1 and stream_id>-1);

    reserve(other_vector.size());

    switch(copy_t)
    {
    case SYNCHRONOUS:
      COPY_FROM<other_device_name, device_name>::execute(other_vector.get_ptr(), this->get_ptr(), current_size);
      break;

    case ASYNCHRONOUS:
      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);

      COPY_FROM<other_device_name, device_name>::execute(other_vector.get_ptr(), this->get_ptr(), current_size, thread_id, stream_id);

      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
  }


  template<typename scalartype, device_type device_name>
  inline typename vector<scalartype, device_name>::vector_scalartype* vector<scalartype, device_name>::get_ptr()
  {
    return data;
  }

  template<typename scalartype, device_type device_name>
  inline typename vector<scalartype, device_name>::vector_scalartype* vector<scalartype, device_name>::get_ptr(int i)
  {
    assert(i>-1 && i<current_size);
    return data+i;
  }

  template<typename scalartype, device_type device_name>
  int& vector<scalartype, device_name>::size(){
    return current_size;
  }

  template<typename scalartype, device_type device_name>
  int& vector<scalartype, device_name>::get_current_size(){
    return current_size;
  }

  template<typename scalartype, device_type device_name>
  int& vector<scalartype, device_name>::get_global_size(){
    return global_size;
  }

  template<typename scalartype, device_type device_name>
  void vector<scalartype, device_name>::reserve(int new_current_size)
  {
    assert(new_current_size>-1);

    if(new_current_size > global_size)
    {
      int new_global_size = (new_current_size/64)*64+64;
      assert(new_global_size >= new_current_size);

      vector_scalartype* new_data = NULL;

      MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_global_size);

      {
        vector_scalartype* tmp_ptr = data;
        data                       = new_data;
        new_data                   = tmp_ptr;
      }

      MEMORY_MANAGEMENT<device_name>::deallocate(new_data);

      global_size  = new_global_size;
      current_size = new_current_size;
    }
    else
      current_size = new_current_size;
  }

  template<typename scalartype, device_type device_name>
  void vector<scalartype, device_name>::resize(int new_current_size)
  {
    assert(new_current_size>-1);

    if(new_current_size > global_size)
    {
      int new_global_size = (new_current_size/64)*64+64;
      assert(new_global_size >= new_current_size);

      vector_scalartype* new_data = NULL;

      MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_global_size);

      COPY_FROM<device_name, device_name>::execute(&data[0], &new_data[0], current_size);

      {
        vector_scalartype* tmp_ptr = data;
        data                       = new_data;
        new_data                   = tmp_ptr;
      }

      MEMORY_MANAGEMENT<device_name>::deallocate(new_data);

      global_size  = new_global_size;
      current_size = new_current_size;
    }
    else
      current_size = new_current_size;
  }

  template<typename scalartype, device_type device_name>
  template<device_type other_device_name>
  scalartype vector<scalartype, device_name>::difference(vector<scalartype, other_device_name>& other_vector)
  {
    if(current_size!=other_vector.size())
      throw std::logic_error(__FUNCTION__);

    vector<scalartype, CPU> cp_this (current_size);
    vector<scalartype, CPU> cp_other(current_size);

    COPY_FROM<      device_name, CPU>::execute(data                  , cp_this .get_ptr(), current_size);
    COPY_FROM<other_device_name, CPU>::execute(other_vector.get_ptr(), cp_other.get_ptr(), current_size);

    scalartype max_dif=0;

    for(int i=0; i<current_size; ++i)
      if(std::fabs(cp_this[i]-cp_other[i]) > max_dif)
        max_dif = std::fabs(cp_this[i]-cp_other[i]);

    if(std::fabs(max_dif)>1.e-6)
      throw std::logic_error(__FUNCTION__);

    return max_dif;
  }

  template<typename scalartype, device_type device_name>
  scalartype vector<scalartype, device_name>::difference(std::vector<scalartype>& other_vector)
  {
    if(current_size!=int(other_vector.size()))
      throw std::logic_error(__FUNCTION__);

    vector<scalartype, CPU> cp_this(current_size);

    COPY_FROM<device_name, CPU>::execute(data, cp_this.get_ptr(), current_size);

    scalartype max_dif=0;

    for(int i=0; i<current_size; ++i)
      if(std::fabs(cp_this[i]-other_vector[i]) > max_dif)
        max_dif = std::fabs(cp_this[i]-other_vector[i]);

    if(std::fabs(max_dif)>1.e-6)
      throw std::logic_error(__FUNCTION__);

    return max_dif;
  }

  template<typename scalartype, device_type device_name>
  void vector<scalartype, device_name>::print()
  {
    MEMORY_MANAGEMENT<device_name>::print(data, current_size, global_size);
  }

}

#endif
