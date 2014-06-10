//-*-C++-*-

#ifndef HDF5_READER_HEADER_H
#define HDF5_READER_HEADER_H

namespace IO
{
  
  /*!
   * \author Peter Staar
   */
  template<>
  class reader<IO::HDF5>
  {
  public:

    typedef H5::H5File file_type;

  public:

    reader();
    ~reader();
    
    bool is_reader() {return true;}
    bool is_writer() {return false;}

    void open_file(std::string file_name);
    void close_file();

    void open_group(std::string name);
    void close_group();

    std::string get_path();

    template<typename arbitrary_struct_t>
    static void from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

    template<typename scalartype>
    void execute(std::string name, scalartype& value);

    template<typename scalar_type>
    void execute(std::string name, std::vector<scalar_type>& value);

    void execute(std::string name, std::string& value);

    void execute(std::string name, std::vector<std::string>& value);

    template<typename domain_type>
    void execute(std::string name, dmn_0<domain_type>& dmn);

    template<typename scalartype, typename domain_type>
    void execute(function<scalartype, domain_type>& f);

    template<typename scalartype, typename domain_type>
    void execute(std::string name, function<scalartype, domain_type>& f);
   
    template<typename scalar_type>
    void execute(std::string name, LIN_ALG::vector<scalar_type , LIN_ALG::CPU>& A);
    
    template<typename scalar_type>
    void execute(std::string name, LIN_ALG::matrix<scalar_type , LIN_ALG::CPU>& A);

  private:

    bool fexists(const char *filename);

  private:

    H5::H5File* my_file;

    std::vector<std::string> my_paths;
  };
  
  reader<IO::HDF5>::reader():
    my_file(NULL),
    my_paths(0)
  {}

  reader<IO::HDF5>::~reader()
  {
    if(my_file != NULL)
      throw std::logic_error(__FUNCTION__);
  }

  void reader<IO::HDF5>::open_file(std::string file_name)
  {
    {// check whether the file exists ...
      std::wifstream tmp(file_name.c_str());
    
      if(!tmp or !tmp.good() or tmp.bad()) 
	{
	  cout << "\n\n\tcannot open file : " << file_name << "\n";
	  throw std::runtime_error(__FUNCTION__);
	}
      else
	{
	  cout << "\n\n\topening file : " << file_name << "\n";
	}
    }

    my_file = new H5File(file_name.c_str(), H5F_ACC_RDONLY);
  }

  void reader<IO::HDF5>::close_file()
  {
    delete my_file;
    my_file = NULL;
  }
  
  void reader<IO::HDF5>::open_group(std::string name)
  {
    my_paths.push_back(name);
  }

  void reader<IO::HDF5>::close_group()
  {
    my_paths.pop_back();
  }
  
  std::string reader<IO::HDF5>::get_path()
  {
    std::string path = "/";

    for(size_t i=0; i<my_paths.size(); i++){
      path = path+my_paths[i];

      if(i<my_paths.size()-1)
	path = path+"/";
    }

    return path;
  }

  template<typename arbitrary_struct_t>
  void reader<IO::HDF5>::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name)
  {
    reader<IO::HDF5> reader_obj;

    reader_obj.open_file(file_name);

    arbitrary_struct.read_write(reader_obj);

    reader_obj.close_file();
  }

  template<typename scalar_type>
  void reader<IO::HDF5>::execute(std::string name, scalar_type& value)//, H5File& file, std::string path)
  {
    std::string full_name = get_path()+"/"+name;
    
    try
      {
	H5::DataSet  dataset   = my_file->openDataSet(full_name.c_str());
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	H5Dread(dataset.getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &value);
      }
    catch(...)
      {
	cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename scalar_type>
  void reader<IO::HDF5>::execute(std::string name, std::vector<scalar_type>& value)//, H5File& file, std::string path)
  {
    std::string full_name = get_path()+"/"+name;
    
    try
      {
	H5::DataSet  dataset   = my_file->openDataSet(full_name.c_str());
    
	value.resize(dataset.getInMemDataSize());
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	H5Dread(dataset.getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &value[0]);
      }
    catch(...)
      {
	cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }

  void reader<IO::HDF5>::execute(std::string name, std::string& value)//, H5File& file, std::string path)
  {
    std::string full_name = get_path()+"/"+name;
    
    try
      {
	H5::DataSet  dataset   = my_file->openDataSet(full_name.c_str());
    
	value.resize(dataset.getInMemDataSize(), 'a');
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	H5Dread(dataset.getId(), IO::HDF5_TYPE<char>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &value[0]);
      }
    catch(...)
      {
	cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }

  void reader<IO::HDF5>::execute(std::string name, std::vector<std::string>& value)//, H5File& file, std::string path)
  {
    try
      {
	open_group(name);

	int size = -1;
	execute("size", size);

	value.resize(size);

	open_group("data");

	for(size_t l=0; l<value.size(); l++)
	  {
	    std::stringstream ss;
	    ss << get_path() << "/" << l;

	    H5::DataSet  dataset   = my_file->openDataSet(ss.str().c_str());
	    
	    value[l].resize(dataset.getInMemDataSize(), 'a');
	    
	    H5::DataSpace dataspace = dataset.getSpace();
	    
	    H5Dread(dataset.getId(), IO::HDF5_TYPE<char>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &value[l][0]);
	  }

	close_group();

	close_group();
      }
    catch(...)
      {
	cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename domain_type>
  void reader<IO::HDF5>::execute(std::string name, dmn_0<domain_type>& dmn)
  {

  }

  template<typename scalartype, typename domain_type>
  void reader<IO::HDF5>::execute(function<scalartype, domain_type>& f)
  {
    execute(f.get_name(), f);
  }

  template<typename scalartype, typename domain_type>
  void reader<IO::HDF5>::execute(std::string name, function<scalartype, domain_type>& f)
  {
    cout << "\n\tstart reading function : " << name;

    try
      {
	open_group(name);

	std::string full_name = get_path()+"/data";

	H5::DataSet  dataset  = my_file->openDataSet(full_name.c_str());
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	H5Dread(dataset.getId(), IO::HDF5_TYPE<scalartype>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &f(0));

	close_group();
      }
    catch(...)
      {
	cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }
  
  template<typename scalar_type>
  void reader<IO::HDF5>::execute(std::string name, LIN_ALG::vector<scalar_type , LIN_ALG::CPU>& V)
  {
    try
      {
	open_group(name);

	std::string full_name = get_path()+"/data";

	H5::DataSet  dataset    = my_file->openDataSet(full_name.c_str());
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	H5Dread(dataset.getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &V[0]);

	close_group();
      }
    catch(...)
      {
	cout << "\n\n\t the vector (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename scalar_type>
  void reader<IO::HDF5>::execute(std::string name, LIN_ALG::matrix<scalar_type , LIN_ALG::CPU>& A)
  {
    try
      {
	open_group(name);

	std::string full_name = get_path()+"/data";

	H5::DataSet  dataset    = my_file->openDataSet(full_name.c_str());
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	H5Dread(dataset.getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &A(0,0));

	close_group();
      }
    catch(...)
      {
	cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }
  

  /*
  template<typename scalar_type>
  void reader<IO::HDF5>::execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& V, H5File& file, std::string path)
  {
    std::string new_path = path+"/"+name;
    
    try
      {
	{
	  int cs=-1;
	  execute("size", cs, file, new_path);
	  
	  V.reserve(cs);
	}

	{
	  std::string  full_path = new_path+"/data";

	  H5::DataSet  dataset    = my_file->openDataSet(full_path.c_str());
	
	  H5::DataSpace dataspace = dataset.getSpace();
	
	  H5Dread(dataset.getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &V[0]);
	}
      }
    catch(...)
      {
	cout << "\n\n\t the vector (" + name + ") does not exist in path : " + path + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename scalar_type>
  void reader<IO::HDF5>::execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A, H5File& file, std::string path)
  {
    std::string new_path = path+"/"+name;
    
    try
      {
	std::vector<int> vec_cs(2,0);
	std::vector<int> vec_gs(2,0);
	
	{
	  execute("current-size", vec_cs, file, new_path);
	  execute("global-size" , vec_gs, file, new_path);
	}

	A.resize_no_copy(std::pair<int,int>(vec_gs[0],vec_gs[1]));
	A.resize_no_copy(std::pair<int,int>(vec_cs[0],vec_cs[1]));

	{
	  std::string  full_path = path+"/"+name+"/data";

	  H5::DataSet  dataset    = my_file->openDataSet(full_path.c_str());
	
	  H5::DataSpace dataspace = dataset.getSpace();
	
	  H5Dread(dataset.getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT, &A(0,0));
	}
      }
    catch(...)
      {
	cout << "\n\n\t the matrix (" + name + ") does not exist in path : " + path + "\n\n";
	//throw std::logic_error(__FUNCTION__);
      }
  }
  */
}

#endif
