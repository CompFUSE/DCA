//-*-C++-*-

#ifndef PACKING_INTERFACE_MPI_H
#define PACKING_INTERFACE_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class packing_interface<MPI_LIBRARY>
  {
  public:

    packing_interface(processor_grouping<MPI_LIBRARY>& grouping_ref);
    ~packing_interface();

    /************************************
     ***  size
     ************************************/

    template<typename scalar_type>
    int get_buffer_size(scalar_type item);

    template<typename scalar_type>
    int get_buffer_size(std::basic_string<scalar_type> str);

    template<typename scalar_type>
    int get_buffer_size(std::vector<scalar_type>& v);

    template<typename scalar_type>
    int get_buffer_size(std::vector<std::basic_string<scalar_type> >& v);

    template<typename scalar_type>
    int get_buffer_size(std::vector<std::vector<scalar_type> >& v);

    template<typename scalar_type, class dmn_type>
    int get_buffer_size(function<scalar_type, dmn_type>& f);

    /************************************
     ***  pack_unpack
     ************************************/

    template<typename object_type>
    void pack_or_unpack(bool packing, int* buffer, int size, int& off_set, object_type& object);

    //private:

    /************************************
     ***  pack
     ************************************/

    template<typename scalar_type>
    void pack(int* buffer, int size, int& off_set, scalar_type& item );

    template<typename scalar_type>
    void pack(int* buffer, int size, int& off_set, std::basic_string<scalar_type>& str );

    template<typename scalar_type>
    void pack(int* buffer, int size, int& off_set, std::vector<scalar_type>& v);

    template<typename scalar_type>
    void pack(int* buffer, int size, int& off_set, std::vector<std::basic_string<scalar_type> >& v);

    template<typename scalar_type>
    void pack(int* buffer, int size, int& off_set, std::vector<std::vector<scalar_type> >& v);

    template<typename scalar_type, class dmn_type>
    void pack(int* buffer, int size, int& off_set, function<scalar_type, dmn_type>& f);

    /************************************
     ***  unpack
     ************************************/

    template<typename scalar_type>
    void unpack(int* buffer, int size, int& off_set, scalar_type& item);

    template<typename scalar_type>
    void unpack(int* buffer, int size, int& off_set, std::basic_string<scalar_type>& str);

    template<typename scalar_type>
    void unpack(int* buffer, int size, int& off_set, std::vector<scalar_type>& v);

    template<typename scalar_type>
    void unpack(int* buffer, int size, int& off_set, std::vector<std::basic_string<scalar_type> >& v);

    template<typename scalar_type>
    void unpack(int* buffer, int size, int& off_set, std::vector<std::vector<scalar_type> >& v);

    template<typename scalar_type, class dmn_type>
    void unpack(int* buffer, int size, int& off_set, function<scalar_type, dmn_type>& f);

  private:

    processor_grouping<MPI_LIBRARY>& grouping;
  };

  packing_interface<MPI_LIBRARY>::packing_interface(processor_grouping<MPI_LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  packing_interface<MPI_LIBRARY>::~packing_interface()
  {}

  template<typename scalar_type>
  int packing_interface<MPI_LIBRARY>::get_buffer_size(scalar_type item)
  {
    int size(0);

    MPI_Pack_size(type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  grouping.get(), &size);

    return size;
  }

  template<typename scalar_type>
  int packing_interface<MPI_LIBRARY>::get_buffer_size(std::basic_string<scalar_type> str)
  {
    /*
      int result = get_buffer_size(str.size());

      int size(0);
      MPI_Pack_size(str.size(), MPI_CHAR, grouping.get(), &size);

      result += size;

      return result;
    */

    int result = get_buffer_size(str.size());

    int count = str.size()*type_map_interface<MPI_LIBRARY, scalar_type>::factor();

    {
      int size(0);
      MPI_Pack_size(count, type_map_interface<MPI_LIBRARY, scalar_type>::value(), grouping.get(), &size);

      result += size;
    }

    return result;
  }

  template<typename scalar_type>
  int packing_interface<MPI_LIBRARY>::get_buffer_size(std::vector<scalar_type>& v)
  {
    int result = get_buffer_size(v.size());

    int count = v.size()*type_map_interface<MPI_LIBRARY, scalar_type>::factor();

    {
      int size(0);
      MPI_Pack_size(count, type_map_interface<MPI_LIBRARY, scalar_type>::value(), grouping.get(), &size);

      result += size;
    }

    return result;
  }

  template<typename scalar_type>
  int packing_interface<MPI_LIBRARY>::get_buffer_size(std::vector<std::basic_string<scalar_type> >& v)
  {
    int result = 0;

    std::vector<int> tmp_sizes(0);
    for(int i=0; i<v.size(); i++)
      tmp_sizes.push_back(v[i].size());

    std::vector<scalar_type> tmp_chars(0);
    for(int i=0; i<v.size(); i++)
      for(int j=0; j<v[i].size(); j++)
        tmp_chars.push_back(v[i][j]);


    result += get_buffer_size(v.size());
    result += get_buffer_size(tmp_sizes);
    result += get_buffer_size(tmp_chars);

    return result;
  }

  template<typename scalar_type>
  int packing_interface<MPI_LIBRARY>::get_buffer_size(std::vector<std::vector<scalar_type> >& v)
  {
    std::vector<scalar_type> tmp;

    tmp.push_back(v.size());
    for(int i=0; i<v.size(); i++)
      {
        tmp.push_back(v[i].size());

        for(int j=0; j<v[i].size(); j++)
          tmp.push_back(v[i][j]);
      }

    return get_buffer_size(tmp);
  }

  template<typename scalar_type, class dmn_type>
  int packing_interface<MPI_LIBRARY>::get_buffer_size(function<scalar_type, dmn_type>& f)
  {
    int result = get_buffer_size(f.size());

    int count = f.size()*type_map_interface<MPI_LIBRARY, scalar_type>::factor();

    {
      int size = 0;
      MPI_Pack_size(count, type_map_interface<MPI_LIBRARY, scalar_type>::value(), grouping.get(), &size);

      result += size;
    }

    return result;
  }

  /************************************
   ***  pack_unpack
   ************************************/

  template<typename object_type>
  void packing_interface<MPI_LIBRARY>::pack_or_unpack(bool packing, int* buffer, int size, int& off_set, object_type& item)
  {
    if(packing)
      pack  (buffer, size, off_set, item);
    else
      unpack(buffer, size, off_set, item);
  }

  /************************************
   ***  pack
   ************************************/

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::pack(int* buffer, int size, int& off_set, scalar_type& item)
  {
    scalar_type* tPtr(&item);

    MPI_Pack(tPtr,
             type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
             type_map_interface<MPI_LIBRARY, scalar_type>::value(),
             buffer,
             size,
             &off_set,
             grouping.get());
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::pack(int* buffer, int size, int& off_set, std::basic_string<scalar_type>& str)
  {
    /*
    // pack the string's length
    int stringSize(str.size());
    pack(buffer,size,off_set,stringSize);

    MPI_Pack(const_cast<char*>(str.c_str()), stringSize, MPI_CHAR,
    buffer, size, &off_set,
    grouping.get());
    */

    // Pack the vector length
    int vectorSize(str.size());
    pack(buffer, size, off_set, vectorSize);

    MPI_Pack(static_cast<scalar_type*>(&str[0]),
             vectorSize*type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
             type_map_interface<MPI_LIBRARY, scalar_type>::value(),
             buffer,
             size,
             &off_set,
             grouping.get());
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::pack(int* buffer, int size, int& off_set, std::vector<scalar_type>& v)
  {
    // Pack the vector length
    int vectorSize(v.size());
    pack(buffer, size, off_set, vectorSize);

    MPI_Pack(static_cast<scalar_type*>(&v[0]),
             vectorSize*type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
             type_map_interface<MPI_LIBRARY, scalar_type>::value(),
             buffer,
             size,
             &off_set,
             grouping.get());
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::pack(int* buffer, int size, int& off_set, std::vector<std::basic_string<scalar_type> >& v)
  {
    std::vector<int> tmp_sizes(0);
    for(int i=0; i<v.size(); i++)
      tmp_sizes.push_back(v[i].size());

    std::vector<scalar_type> tmp_chars(0);
    for(int i=0; i<v.size(); i++)
      for(int j=0; j<v[i].size(); j++)
        tmp_chars.push_back(v[i][j]);

    {
      // Pack the vector length
      int v_size = v.size();
      pack(buffer, size, off_set, v_size);

      pack(buffer, size, off_set, tmp_sizes);
      pack(buffer, size, off_set, tmp_chars);
    }
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::pack(int* buffer, int size, int& off_set, std::vector<std::vector<scalar_type> >& v)
  {
    std::vector<scalar_type> tmp;

    tmp.push_back(v.size());
    for(int i=0; i<v.size(); i++)
      {
        tmp.push_back(v[i].size());

        for(int j=0; j<v[i].size(); j++)
          tmp.push_back(v[i][j]);
      }

    pack(buffer, size, off_set, tmp);
  }

  template<typename scalar_type, class dmn_type>
  void packing_interface<MPI_LIBRARY>::pack(int* buffer, int size, int& off_set, function<scalar_type, dmn_type>& f)
  {
    // Pack the vector length
    int function_size(f.size());
    pack(buffer, size, off_set, function_size);

    MPI_Pack(static_cast<scalar_type*>(&f(0)),
             function_size*type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
             type_map_interface<MPI_LIBRARY, scalar_type>::value(),
             buffer,
             size,
             &off_set,
             grouping.get());

  }

  /************************************
   ***  unpack
   ************************************/

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::unpack(int* buffer, int size, int& off_set, scalar_type& item)
  {
    scalar_type tmp;

    MPI_Unpack(buffer,
               size,
               &off_set,
               &tmp,
               type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
               type_map_interface<MPI_LIBRARY, scalar_type>::value(),
               grouping.get());

    item = tmp;
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::unpack(int* buffer, int size, int& off_set, std::basic_string<scalar_type>& str)
  {
    /*
    // Unpack the string length
    int stringSize(0);
    unpack(buffer, size, off_set, stringSize);

    char stringBuffer[stringSize];
    MPI_Unpack(buffer,
    size,
    &off_set,
    stringBuffer,
    stringSize,
    MPI_CHAR,
    grouping.get());

    std::basic_string<scalar_type> tmp(stringBuffer,stringSize);

    str = tmp;
    */

    // UnPack the vector length
    int vectorSize(0);
    unpack(buffer, size, off_set, vectorSize);

    str.resize(vectorSize);

    // UnPack the vector
    MPI_Unpack(buffer,
               size,
               &off_set,
               static_cast<scalar_type*>(&str[0]),
               type_map_interface<MPI_LIBRARY, scalar_type>::factor()*vectorSize,
               type_map_interface<MPI_LIBRARY, scalar_type>::value(),
               grouping.get());
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::unpack(int* buffer, int size, int& off_set, std::vector<scalar_type>& v)
  {
    // UnPack the vector length
    int vectorSize(0);
    unpack(buffer, size, off_set, vectorSize);

    v.resize(vectorSize);

    // UnPack the vector
    MPI_Unpack(buffer,
               size,
               &off_set,
               static_cast<scalar_type*>(&v[0]),
               type_map_interface<MPI_LIBRARY, scalar_type>::factor()*vectorSize,
               type_map_interface<MPI_LIBRARY, scalar_type>::value(),
               grouping.get());
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::unpack(int* buffer, int size, int& off_set, std::vector<std::basic_string<scalar_type> >& v)
  {
    int v_size;

    std::vector<int >        tmp_sizes(0);
    std::vector<scalar_type> tmp_chars(0);

    {
      unpack(buffer, size, off_set, v_size);

      unpack(buffer, size, off_set, tmp_sizes);
      unpack(buffer, size, off_set, tmp_chars);
    }

    {
      v.resize(v_size);

      for(int i=0; i<v.size(); i++)
        v[i].resize(tmp_sizes[i]);

      int index=0;
      for(int i=0; i<v.size(); i++){
        for(int j=0; j<v[i].size(); j++){
          v[i][j] = tmp_chars[index];

          index+=1;
        }
      }
    }
  }

  template<typename scalar_type>
  void packing_interface<MPI_LIBRARY>::unpack(int* buffer, int size, int& off_set, std::vector<std::vector<scalar_type> >& v)
  {
    std::vector<scalar_type> tmp;

    unpack(buffer, size, off_set, tmp);

    v.resize(0);

    for(int i=1; i<tmp.size(); ){

      int v_size = tmp[i];

      std::vector<scalar_type> tmp_i(0);
      for(int j=i+1; j<i+1+v_size; j++)
        tmp_i.push_back(tmp[j]);

      v.push_back(tmp_i);

      i += (v_size+1);
    }
  }

  template<typename scalar_type, class dmn_type>
  void packing_interface<MPI_LIBRARY>::unpack(int* buffer, int size, int& off_set, function<scalar_type, dmn_type>& f)
  {
    // UnPack the vector length
    int function_size(0);
    unpack(buffer, size, off_set, function_size);

    // UnPack the vector
    MPI_Unpack(buffer,
               size,
               &off_set,
               static_cast<scalar_type*>(&f(0)),
               type_map_interface<MPI_LIBRARY, scalar_type>::factor()*function_size,
               type_map_interface<MPI_LIBRARY, scalar_type>::value(),
               grouping.get());
  }

}

#endif
