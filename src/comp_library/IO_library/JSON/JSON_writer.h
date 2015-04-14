//-*-C++-*-

#ifndef JSON_WRITER_HEADER_H
#define JSON_WRITER_HEADER_H

namespace IO
{

  /*!
   *  \author Peter Staar
   */
  template<>
  class writer<IO::JSON>
  {
  public:

    typedef std::stringstream file_type;

    typedef JSONPARSER::Whatever     JsonAccessor;
    typedef JSONPARSER::JSON_context JsonDataType;

  public:

    writer();
    ~writer();

    bool is_reader() {return false;}
    bool is_writer() {return true;}

    file_type& open_file(std::string file_name_ref, bool overwrite=true);
    void       close_file();

    void open_group(std::string new_path);
    void close_group();

    std::string get_path();

    template<typename arbitrary_struct_t>
    static void to_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

    template<typename scalartype>
    void execute(std::string name,             scalartype&  value);

    template<typename s_t_0, typename s_t_1>
    void execute(std::string name, std::pair<s_t_0, s_t_1>& value);

    template<typename scalartype>
    void execute(std::string name, std::vector<scalartype>& value);

    template<typename scalartype>
    void execute(std::string name, std::vector<std::vector<scalartype> >& value);

    void execute(std::string name,             std::string&  value);

    void execute(std::string name, std::vector<std::string>& value);

    template<typename domain_type>
    void execute(std::string name,  dmn_0<domain_type>& dmn);

    template<typename scalar_type, typename domain_type>
    void execute(FUNC_LIB::function<scalar_type, domain_type>& f);

    template<typename scalar_type, typename domain_type>
    void execute(std::string name, FUNC_LIB::function<scalar_type, domain_type>& f);

    template<typename scalar_type, typename domain_type>
    void execute(FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f);

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

    template<class stream_type>
    static void execute(stream_type& ss, const JsonAccessor& parseResult);

  private:

    std::stringstream ss;

    std::string file_name;

    std::string path;

    std::vector<int> elements_in_group;
  };

  writer<IO::JSON>::writer():
    file_name(""),
    path     (""),
    elements_in_group(0)
  {
    ss << std::fixed;
    ss.precision(16);
  }

  writer<IO::JSON>::~writer()
  {}

  std::stringstream& writer<IO::JSON>::open_file(std::string my_file_name, bool overwrite)
  {
    file_name = my_file_name;

    ss << "{\n";

    elements_in_group.push_back(0);

    return ss;
  }

  void writer<IO::JSON>::close_file()
  {
    ss << "\n}";

    {
      std::ofstream of(&file_name[0]);

      of << ss.str();

      of.flush();
      of.close();
    }
  }

  void writer<IO::JSON>::open_group(std::string name)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n\n";

    ss << get_path() << "\"" << name << "\"" << " : \n";

    ss << get_path() << "{\n";

    elements_in_group.push_back(0);
  }

  void writer<IO::JSON>::close_group()
  {
    elements_in_group.pop_back();

    ss << "\n" << get_path() << "}";

    elements_in_group.back() += 1;
  }

  std::string writer<IO::JSON>::get_path()
  {
    std::stringstream ss;
    for(size_t i=0; i<elements_in_group.size(); i++)
      ss << "\t";

    return ss.str();
  }

  template<typename arbitrary_struct_t>
  void writer<IO::JSON>::to_file(arbitrary_struct_t& arbitrary_struct, std::string file_name)
  {
    writer<IO::JSON> wr_obj;

    wr_obj.open_file(file_name);

    arbitrary_struct.read_write(wr_obj);

    wr_obj.close_file();
  }

  template<typename scalartype>
  void writer<JSON>::execute(std::string name, scalartype& value)//, file_type& ss), std::string path, bool is_ending)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n";

    ss << get_path() << "\"" << name << "\" : " << value;

    elements_in_group.back() += 1;
  }

  template<typename s_t_0, typename s_t_1>
  void writer<JSON>::execute(std::string name, std::pair<s_t_0, s_t_1>& value)//, file_type& ss), std::string path, bool is_ending)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n";

    ss << get_path() << "\"" << name << "\" : [" << value.first << ", " << value.second << "]";

    elements_in_group.back() += 1;
  }

  template<typename scalartype>
  void writer<JSON>::execute(std::string name, std::vector<scalartype>& value)//, file_type& ss)//, std::string path, bool is_ending)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n";

    ss << get_path() << "\"" << name << "\" : [";

    for(size_t i=0; i<value.size(); i++)
      {
        ss << value[i];

        if(i < value.size()-1)
          ss << ", ";
      }

    ss << "]";

    elements_in_group.back() += 1;
  }

  template<typename scalartype>
  void writer<JSON>::execute(std::string name, std::vector<std::vector<scalartype> >& value)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n";

    ss << get_path() << "\"" << name << "\" : [";

    std::string indent="";
    indent.resize(name.size()+6, ' ');

    for(size_t i=0; i<value.size(); i++)
      {
        ss << "[";

        for(size_t j=0; j<value[i].size(); j++)
          {
            ss << value[i][j];

            if(j < value[i].size()-1)
              ss << ", ";
          }
        ss << "]";

        if(i < value.size()-1)
          ss << ",\n" << get_path() << indent ;
      }
    ss << "]";

    elements_in_group.back() += 1;
  }

  void writer<JSON>::execute(std::string name, std::string& value)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n";

    ss << get_path() << "\"" << name << "\" : \"" << value << "\"";

    elements_in_group.back() += 1;
  }

  void writer<JSON>::execute(std::string name, std::vector<std::string>& value)
  {
    if(elements_in_group.back() != 0)
      ss << ",\n";

    ss << get_path() << "\"" << name << "\" : [";

    for(size_t i=0; i<value.size(); i++)
      {
        ss << "\"" << value[i] << "\"";

        if(i == value.size()-1)
          ss << "]";
        else
          ss << ", ";
      }

    elements_in_group.back() += 1;
  }

  template<typename domain_type>
  void writer<IO::JSON>::execute(std::string name, dmn_0<domain_type>& dmn)
  {
    open_group(name);

    execute("name", dmn.get_name());

    execute("elements", dmn.get_elements());

    close_group();

    elements_in_group.back() += 1;
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::JSON>::execute(FUNC_LIB::function<scalar_type, domain_type>& f)
  {
    cout << "\t starts writing function : " << f.get_name() << "\n";

    execute(f.get_name(), f);
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::JSON>::execute(std::string name, FUNC_LIB::function<scalar_type, domain_type>& f)
  {
    open_group(name);

    execute("name", f.get_name());

    {
      std::vector<int> vec(0);

      for(int l=0; l<f.signature(); l++)
        vec.push_back(f[l]);

      execute("domain-sizes", vec);//, file, new_path);
    }

    ss << ",\n\n" << get_path() << "\"data\" : [";

//     ss << std::fixed;
//     ss.precision(16);

    int* subind = new int[f.signature()];
    for(int i=0; i<f.size(); i++)
      {
        ss << "[";

        int dummy = i;
        f.linind_2_subind(dummy,subind);
        for(int j=0; j<f.signature(); j++)
          {
            ss << (subind[j]);
            ss<< ", ";
          }

        ss << f(i);

        if(i == f.size()-1)
          ss << "]";
        else
          ss << "],\n" << get_path() << "          ";
      }
    delete [] subind;

    ss << "]";

    close_group();

    elements_in_group.back() += 1;
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::JSON>::execute(FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f)
  {
    cout << "\t starts writing function : " << f.get_name() << "\n";

    execute(f.get_name(), f);
  }

  template<typename scalar_type, typename domain_type>
  void writer<IO::JSON>::execute(std::string name, FUNC_LIB::function<std::complex<scalar_type>, domain_type>& f)
  {
    open_group(name);

    execute("name", f.get_name());

    {
      std::vector<int> vec(0);

      for(int l=0; l<f.signature(); l++)
        vec.push_back(f[l]);

      execute("domain-sizes", vec);
    }

    ss << ",\n\n" << get_path() << "\"data\" : [";

//     ss << std::fixed;
//     ss.precision(16);

    int* subind = new int[f.signature()];
    for(int i=0; i<f.size(); i++)
      {
        ss << "[";

        int dummy = i;
        f.linind_2_subind(dummy,subind);
        for(int j=0; j<f.signature(); j++)
          {
            ss << (subind[j]);
            ss<< ", ";
          }

        ss << real(f(i)) << ", " << imag(f(i));

        if(i == f.size()-1)
          ss << "]";
        else
          ss << "],\n" << get_path() << "          ";
      }
    delete [] subind;

    ss << "]";

    close_group();

    elements_in_group.back() += 1;
  }

  template<typename scalar_type>
  void writer<IO::JSON>::execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& V)
  {
    open_group(name);

    execute("name", V.get_name());

    execute("current-size", V.get_current_size());
    execute("global-size" , V.get_global_size());

//     ss << std::fixed;
//     ss.precision(16);

    {
      ss << ",\n\n" << get_path() << "\"data\" : [";

      for(int j=0; j<V.get_current_size()-1; j++)
        ss << V[j] << ", ";

      ss << V[V.get_current_size()-1] << "]";
    }

    close_group();

    elements_in_group.back() += 1;
  }

  template<typename scalar_type>
  void writer<IO::JSON>::execute(std::string name, LIN_ALG::vector<std::complex<scalar_type>, LIN_ALG::CPU>& V)
  {
    open_group(name);

    execute("name", V.get_name());

    execute("current-size", V.get_current_size());
    execute("global-size" , V.get_global_size());

//     ss << std::fixed;
//     ss.precision(16);

    {
      ss << ",\n\n" << get_path() << "\"data-real\" : [";

      for(int j=0; j<V.get_current_size()-1; j++)
        ss << real(V[j]) << ", ";

      ss << real(V[V.get_current_size()-1]) << "]";
    }

    {
      ss << ",\n\n" << get_path() << "\"data-imag\" : [";

      for(int j=0; j<V.get_current_size()-1; j++)
        ss << imag(V[j]) << ", ";

      ss << imag(V[V.get_current_size()-1]) << "]";
    }

    close_group();

    elements_in_group.back() += 1;
  }

  template<typename scalar_type>
  void writer<IO::JSON>::execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A)
  {
    open_group(name);

    execute("name", A.get_name());

    execute("current-size", A.get_current_size());
    execute("global-size" , A.get_global_size());

    ss << ",\n\n" << get_path() << "\"data\" : [";

//     ss << std::fixed;
//     ss.precision(16);

    for(int i=0; i<A.get_current_size().first; i++)
      {
        ss << "[";

        for(int j=0; j<A.get_current_size().second-1; j++)
          ss << A(i,j) << ", ";

        if(i == A.get_current_size().first-1)
          ss << A(i, A.get_current_size().second-1) << "]";
        else
          ss << A(i, A.get_current_size().second-1) << "],\n" << get_path() << "          ";
      }

    ss << "]";

    close_group();

    elements_in_group.back() += 1;
  }

  template<typename scalar_type>
  void writer<IO::JSON>::execute(std::string name, LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>& A)
  {
    open_group(name);

    execute("name", A.get_name());

    execute("current-size", A.get_current_size());
    execute("current-size", A.get_global_size());

    {
      ss << ",\n\n" << get_path() << "\"data-real\" : [";

//       ss << std::fixed;
//       ss.precision(16);

      for(int i=0; i<A.get_current_size().first; i++)
        {
          ss << "[";

          for(int j=0; j<A.get_current_size().second-1; j++)
            ss << real(A(i,j)) << ", ";

          if(i == A.get_current_size().first-1)
            ss << real(A(i, A.get_current_size().second-1)) << "]";
          else
            ss << real(A(i, A.get_current_size().second-1)) << "],\n" << get_path() << "               ";
        }

      ss << "]";
    }

    {
      ss << ",\n\n" << get_path() << "\"data-imag\" : [";

//       ss << std::fixed;
//       ss.precision(16);

      for(int i=0; i<A.get_current_size().first; i++)
        {
          ss << "[";

          for(int j=0; j<A.get_current_size().second-1; j++)
            ss << imag(A(i,j)) << ", ";

          if(i == A.get_current_size().first-1)
            ss << imag(A(i, A.get_current_size().second-1)) << "]";
          else
            ss << imag(A(i, A.get_current_size().second-1)) << "],\n" << get_path() << "               ";
        }

      ss << "]";
    }

    close_group();

    elements_in_group.back() += 1;
  }










  template<class stream_type>
  void writer<IO::JSON>::execute(stream_type& os, const JsonAccessor& w)
  {
    static int level = -1;

    //typedef          std::map<std::wstring, JSONPARSER::Whatever>                 WhateverMap;
    typedef typename std::map<std::wstring, JSONPARSER::Whatever>::const_iterator WhateverMapItr;

    //typedef std::string                  StringType;

    switch (w.type)
      {
      case JSONPARSER::WHATEVER_MAT:
        {
          std::string wfilename(w.filename.begin(),w.filename.end());
          os << "{ 'fileName': '" << wfilename << "'"
             << ", 'startPos': "  << w.startPos
             << ", 'endPos': "    << w.endPos << "}";
          break;
        }

      case JSONPARSER::WHATEVER_MAP:
        {
          level += 1;

          os << "\n";
          for(int l=0; l<level; l++)
            os << "\t";
          os << "{\n";

          int index = 0;
          for(WhateverMapItr itr = w.whateverMap.begin(); itr != w.whateverMap.end(); itr++)
            {
              const std::wstring& wkey = (*itr).first;
              const std::string    key(wkey.begin(), wkey.end());

              for(int l=0; l<level; l++)
                os << "\t";

              os << "\"" << key << "\" : ";// << (*itr).second;

              writer<IO::JSON>::execute(os, (*itr).second);

              if(int(w.whateverMap.size()) == index+1)
                os << "";
              else
                os << ",\n";

              index += 1;
            }

          os << "\n";
          for(int l=0; l<level; l++)
            os << "\t";
          os << "}";

          level -= 1;

          break;
        }

      case JSONPARSER::WHATEVER_VECTOR:
        {
          os << "[";
          for(size_t i=0; i<w.whateverVector.size(); i++)
            {
              writer<IO::JSON>::execute(os, w.whateverVector[i]);

              //os << w.whateverVector[i];
              if(i<w.whateverVector.size()-1)
                os << ", ";
            }

          os << "]";

          break;
        }

      case JSONPARSER::WHATEVER_MATRIX:
        os << "WHATEVER_MATRIX";
        break;

      case JSONPARSER::WHATEVER_STRING:
        {
          const std::string tmp(w.valueString.begin(), w.valueString.end());
          os << "\"" << tmp << "\"";
        }
        break;

      case JSONPARSER::WHATEVER_INTEGER:
        os << w.whateverInteger;
        break;

      case JSONPARSER::WHATEVER_DOUBLE:
        os << w.whateverDouble;
        break;

      case JSONPARSER::WHATEVER_BOOL:
        os << w.whateverBool;
        break;

      case JSONPARSER::WHATEVER_UNKNOWN:
        os <<"WHATEVER_UNKNOWN";
        break;

      default:
        throw std::logic_error("typeName given wrong type");
      }
  }
}

#endif
