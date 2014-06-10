//-*-C++-*-

#ifndef HDF5_GROUP_HEADER_H
#define HDF5_GROUP_HEADER_H

namespace IO
{
  class hdf5_group
  {
  public:
    
    hdf5_group(hid_t file_id, const std::string path);
    ~hdf5_group();
    
    hid_t get_id();

  private:
    
    hid_t id;    
  };
  
  hdf5_group::hdf5_group(hid_t file_id, const std::string path)
  {
    id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
   
    if(id < 0) 
      throw std::logic_error(__FUNCTION__);
  }
  
  hdf5_group::~hdf5_group()
  {            
    if(H5Gclose(id) < 0) 
      throw std::logic_error(__FUNCTION__);
  }
 
  hid_t hdf5_group::get_id()
  {
    return id;
  }
 
}

#endif
