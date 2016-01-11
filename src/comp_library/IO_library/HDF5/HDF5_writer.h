//-*-C++-*-

#ifndef HDF5_WRITER_HEADER_H
#define HDF5_WRITER_HEADER_H

#include <vector>

namespace IO
{

  /*!
   * \author Peter Staar
   */
  template<>
  class writer<IO::HDF5>
  {
  public:

    typedef H5::H5File file_type;

  public:

    writer();
    ~writer();

    bool is_reader() {return false;}
    bool is_writer() {return true;}

    H5::H5File& open_file(std::string file_name_ref, bool overwrite=true);
    void close_file();

    void open_group(std::string new_path);
    void close_group();

    std::string get_path();

    template<typename arbitrary_struct_t>
    static void to_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

    template<typename scalar_type>
    void execute(std::string name,             scalar_type   value);

    template<typename scalar_type>
    void execute(std::string name, std::pair<scalar_type, scalar_type>& value);

    template<typename scalar_type>
    void execute(std::string name, std::vector<scalar_type>& value);

    template<typename scalar_type>
    void execute(std::string name, std::vector< std::complex<scalar_type> >& value);

    void execute(std::string name,             std::string   value);

    void execute(std::string name, std::vector<std::string>& value);

    template<typename scalar_type>
    void execute(std::string name, std::vector<std::vector<scalar_type> >& value);

    template<typename domain_type>
    void execute(std::string name, dmn_0<domain_type>& dmn);

    template<typename scalar_type, typename domain_type>
    void execute(FUNC_LIB::function<             scalar_type , domain_type>& f);

    template<typename scalar_type, typename domain_type>
    void execute(FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f);

    template<typename scalar_type, typename domain_type>
    void execute(std::string name, FUNC_LIB::function<             scalar_type , domain_type>& f);

    template<typename scalar_type, typename domain_type>
    void execute(std::string name, FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f);

    template<typename scalar_type>
    void execute(std::string name, LIN_ALG::vector<             scalar_type , LIN_ALG::CPU>& A);

    template<typename scalar_type>
    void execute(std::string name, LIN_ALG::vector<std::complex<scalar_type>, LIN_ALG::CPU>& A);

    template<typename scalar_type>
    void execute(std::string name, LIN_ALG::matrix<             scalar_type , LIN_ALG::CPU>& A);

    template<typename scalar_type>
    void execute(std::string name, LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>& A);

  private:

    bool fexists(const char *filename);

  private:

    H5::H5File* my_file;

    hid_t       file_id;

    std::vector<H5::Group* >  my_group;
    std::vector<std::string>  my_paths;
  };

  writer<IO::HDF5>::writer():
    my_file(NULL),
    my_group(0),
    my_paths(0)
  {}

  writer<IO::HDF5>::~writer()
  {
    if(my_file != NULL or my_group.size() != 0)
      throw std::logic_error(__FUNCTION__);
  }

  bool writer<IO::HDF5>::fexists(const char *filename)
  {
    std::ifstream ifile(filename);
    return bool(ifile);
  }

  H5::H5File& writer<IO::HDF5>::open_file(std::string file_name, bool overwrite)
  {
    if(my_file != NULL or my_group.size() != 0)
      throw std::logic_error(__FUNCTION__);

    if(overwrite)
      {
	file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      }
    else
      {
	if(fexists(file_name.c_str()))
	  file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	else
	  file_id = H5Fcreate(file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
      }

    my_file = new H5File(file_name.c_str(), H5F_ACC_RDWR);

    return (*my_file);
  }

  void writer<IO::HDF5>::close_file()
  {
    delete my_file;
    my_file = NULL;

    H5Fclose(file_id);
  }

  void writer<IO::HDF5>::open_group(std::string name)
  {
    my_paths.push_back(name);
    my_group.push_back(NULL);

    my_group.back() = new Group(my_file->createGroup(get_path().c_str()));
  }

  void writer<IO::HDF5>::close_group()
  {
    delete my_group.back();

    my_group.pop_back();
    my_paths.pop_back();
  }

  std::string writer<IO::HDF5>::get_path()
  {
    std::string path = "/";

    for(size_t i=0; i<my_paths.size(); i++){
      path = path+my_paths[i];

      if(i<my_paths.size()-1)
        path = path+"/";
    }

    //cout << path << endl;

    return path;
  }

  template<typename arbitrary_struct_t>
  void writer<IO::HDF5>::to_file(arbitrary_struct_t& arbitrary_struct, std::string file_name)
  {
    writer<IO::HDF5> wr_obj;

    wr_obj.open_file(file_name);

    arbitrary_struct.read_write(wr_obj);

    wr_obj.close_file();
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, scalar_type value)
  {
    H5File& file     = (*my_file);
    std::string path = get_path();

    hsize_t  dims[1];

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dims[0] = 1;
      dataspace = new DataSpace(1, dims);

      std::string full_name = path+"/"+name;
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &value);
    }

    delete dataset;
    delete dataspace;
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, std::pair<scalar_type, scalar_type>& value)
  {
    H5File& file     = (*my_file);
    std::string path = get_path();

    hsize_t  dims[1];

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dims[0] = 2;
      dataspace = new DataSpace(1, dims);

      std::string full_name = path+"/"+name;
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &(value.first));
    }

    delete dataset;
    delete dataspace;
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, std::vector<scalar_type>& value)//, H5File& file, std::string path)
  {
    if(value.size()>0)
      {
        H5File& file     = (*my_file);
        std::string path = get_path();

        hsize_t  dims[1];

        DataSet*   dataset   = NULL;
        DataSpace* dataspace = NULL;

        {
          dims[0] = value.size();
          dataspace = new DataSpace(1, dims);

          std::string full_name = path+"/"+name;
          dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

          H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &value[0]);
        }

        delete dataset;
        delete dataspace;
      }
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, std::vector< std::complex<scalar_type> >& value)
  {
    H5File& file     = (*my_file);
    std::string path = get_path();

    hsize_t  dims[2];

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dims[0] = 2;
      dims[1] = value.size();
      dataspace = new DataSpace(2, dims);

      std::string full_name = path+"/"+name;
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &value[0]);
    }

    delete dataset;
    delete dataspace;
  }

  void writer<IO::HDF5>::execute(std::string name, std::string value)//, H5File& file, std::string path)
  {
    if(value.size()>0)
      {
        H5File& file     = (*my_file);
        std::string path = get_path();

        hsize_t  dims[1];

        DataSet*   dataset   = NULL;
        DataSpace* dataspace = NULL;

        {
          dims[0] = value.size();
          dataspace = new DataSpace(1, dims);

          std::string full_name = path+"/"+name;
          dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<char>::get_PredType(), *dataspace));

          H5Dwrite(dataset->getId(), IO::HDF5_TYPE<char>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &value[0]);
        }

        delete dataset;
        delete dataspace;
      }
  }

  void writer<IO::HDF5>::execute(std::string name, std::vector<std::string>& value)//, H5File& file, std::string path)
  {
    if(value.size()>0)
      {
        H5File& file     = (*my_file);

        open_group(name);

        execute("size", value.size());//, file, new_path);

        open_group("data");

        hsize_t  dims[1];

        DataSet*   dataset   = NULL;
        DataSpace* dataspace = NULL;

        for(size_t l=0; l<value.size(); l++)
          {
            dims[0] = value[l].size();
            dataspace = new DataSpace(1, dims);

            std::stringstream ss;
            ss << get_path() << "/" <<l;

            dataset = new DataSet(file.createDataSet(ss.str().c_str(), IO::HDF5_TYPE<char>::get_PredType(), *dataspace));

            H5Dwrite(dataset->getId(), IO::HDF5_TYPE<char>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &(value[l][0]));
          }

        close_group();

        delete dataset;
        delete dataspace;

        close_group();
      }
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, std::vector<std::vector<scalar_type> >& value)
  {
    if(value.size()>0)
      {

        H5File& file = (*my_file);

        bool all_the_same_size=true;

        std::vector<size_t> dim(2,0);
        {
          dim[0] = value   .size();
          dim[1] = value[0].size();

          for(size_t l=0; l<dim[0]; l++)
            if(dim[1] != value[l].size())
              all_the_same_size = false;
        }

        open_group(name);

        execute("equal-size", all_the_same_size);

        DataSet*   dataset   = NULL;
        DataSpace* dataspace = NULL;

        if(all_the_same_size)
          {
            execute("size"      , dim);

            hsize_t  dims[2];

            {
              dims[0] = value   .size();
              dims[1] = value[0].size();

              dataspace = new DataSpace(2, dims);

              std::string full_name = get_path()+"/data";
              dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

              scalar_type* tmp = new scalar_type[dims[0]*dims[1]];

              /*
                for(int i=0; i<dims[0]; i++)
                for(int j=0; j<dims[1]; j++)
                tmp[i+j*dims[0]] = value[i][j];
              */

              // hdf5 has row-major ordering!
              for(int i=0; i<dims[0]; i++)
                for(int j=0; j<dims[1]; j++)
                  tmp[i*dims[1]+j] = value[i][j];

              H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, tmp);

              delete [] tmp;
            }
          }
        else
          {
            execute("size_i", dim[0]);

            dim.resize(0);
            for(size_t l=0; l<value.size(); l++)
              dim.push_back(value[l].size());

            execute("size_j", dim);

            hsize_t  dims[1];

            open_group("data");
            for(size_t l=0; l<value.size(); l++)
              {
                dims[0] = value[l].size();
                dataspace = new DataSpace(1, dims);

                std::stringstream ss;
                ss << get_path() << "/" << l;

                dataset   = new DataSet(file.createDataSet(ss.str().c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

                H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &value[l]);
              }
            close_group();
          }

        delete dataset;
        delete dataspace;

        close_group();
      }
  }

  template<typename domain_type>
  void writer<IO::HDF5>::execute(std::string name, dmn_0<domain_type>& dmn)
  {
    open_group(name);

    execute("name"    , dmn.get_name());
    execute("elements", dmn.get_elements());

    close_group();
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::HDF5>::execute(FUNC_LIB::function<scalar_type, domain_type>& f)
  {
    if(f.size()==0)
      return;

    std::cout << "\t starts writing function : " << f.get_name() << "\n";

    execute(f.get_name(), f);
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::HDF5>::execute(FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f)
  {
    if(f.size()==0)
      return;

    std::cout << "\t starts writing function : " << f.get_name() << "\n";

    execute(f.get_name(), f);
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::HDF5>::execute(std::string name, FUNC_LIB::function<scalar_type, domain_type>& f)
  {
    if(f.size()==0)
      return;

    H5File& file = (*my_file);

    open_group(name);

    std::string new_path = get_path();

    {
      execute("name", f.get_name());

      {
        std::vector<int> vec(0);

        for(int l=0; l<f.signature(); l++)
          vec.push_back(f[l]);

        execute("domain-sizes", vec);
      }

      {
        int N_dmns = f.signature();

        std::vector<hsize_t> dims(N_dmns);  // hsize_t dims[N_dmns];

        //      for(int l=0; l<N_dmns; l++)
        //        dims[l] = f[l];

        // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
        for(int l=0; l<N_dmns; l++)
          dims[N_dmns-1-l] = f[l];

        DataSet*   dataset   = NULL;
        DataSpace* dataspace = NULL;

        {
          dataspace = new DataSpace(N_dmns, &dims[0]);

          std::string full_name = new_path+"/data";
          dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

          H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &f(0));
        }

        delete dataset;
        delete dataspace;
      }
    }

    close_group();
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::HDF5>::execute(std::string name, FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f)
  {
    if(f.size()==0)
      return;

    H5File& file = (*my_file);

    open_group(name);

    {
      execute("name", f.get_name());

      {
        std::vector<int> vec(0);

        for(int l=0; l<f.signature(); l++)
          vec.push_back(f[l]);

        execute("domain-sizes", vec);
      }

      {
        int N_dmns = f.signature();

        std::vector<hsize_t> dims(N_dmns+1);  // hsize_t dims[N_dmns+1];

        //      dims[0] = 2;
        //      for(int l=0; l<N_dmns; l++)
        //        dims[1+l] = f[l];

        // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
        dims[N_dmns] = 2;
        for(int l=0; l<N_dmns; l++)
          dims[N_dmns-1-l] = f[l];

        DataSet*   dataset   = NULL;
        DataSpace* dataspace = NULL;

        {
          dataspace = new DataSpace(N_dmns+1, &dims[0]);

          std::string full_name = get_path()+"/data";
          dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

          H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &f(0));
        }

        delete dataset;
        delete dataspace;
      }
    }

    close_group();
  }


  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& V)
  {
    H5File& file = (*my_file);

    open_group(name);

    execute("name", V.get_name());
    execute("size", V.get_size());

    hsize_t dims[1];

    dims[0] = V.size();

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dataspace = new DataSpace(1, dims);

      std::string full_name = get_path()+"/data";
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &V[0]);
    }

    delete dataset;
    delete dataspace;

    close_group();
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, LIN_ALG::vector<std::complex<scalar_type>, LIN_ALG::CPU>& V)
  {
    H5File& file = (*my_file);

    open_group(name);

    execute("name", V.get_name());
    execute("size", V.get_current_size());

    hsize_t dims[2];

    dims[0] = 2;
    dims[1] = V.get_current_size();

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dataspace = new DataSpace(2, dims);

      std::string full_name = get_path()+"/data";
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &V[0]);
    }

    delete dataset;
    delete dataspace;

    close_group();
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A)
  {
    H5File& file = (*my_file);

    open_group(name);

    {
      execute("name", A.get_name());

      execute("current-size", A.get_current_size());
      execute("global-size" , A.get_global_size());
    }

    hsize_t dims[2];

    dims[0] = A.get_global_size().first;
    dims[1] = A.get_global_size().second;

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dataspace = new DataSpace(2, dims);

      std::string full_name = get_path()+"/data";
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &A(0,0));
    }

    delete dataset;
    delete dataspace;

    close_group();
  }

  template<typename scalar_type>
  void writer<IO::HDF5>::execute(std::string name, LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>& A)
  {
    H5File& file = (*my_file);

    open_group(name);

    {
      execute("name", A.get_name());

      std::vector<int> vec(2,0);

      vec[0] = A.get_current_size().first;
      vec[1] = A.get_current_size().second;
      execute("current-size", vec);

      vec[0] = A.get_global_size().first;
      vec[1] = A.get_global_size().second;
      execute("global-size" , vec);
    }

    hsize_t dims[3];

    dims[0] = 2;
    dims[1] = A.get_global_size().first;
    dims[2] = A.get_global_size().second;

    DataSet*   dataset   = NULL;
    DataSpace* dataspace = NULL;

    {
      dataspace = new DataSpace(3, dims);

      std::string full_name = get_path()+"/data";
      dataset   = new DataSet(file.createDataSet(full_name.c_str(), IO::HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), IO::HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT, &A(0,0));
    }

    delete dataset;
    delete dataspace;

    close_group();
  }

}

#endif
