//-*-C++-*-

#ifndef HDF5_DATASET_HEADER_H
#define HDF5_DATASET_HEADER_H

namespace IO
{
  class hdf5_dataset
  {
  public:
        
    hdf5_dataset(hid_t group_id, const std::string& name);
    
    ~hdf5_dataset();
    
    hid_t get_id();

  private:
    
    hid_t id;
  };
  
  hdf5_dataset::hdf5_dataset(hid_t group_id, const std::string& name)
  {
    id = H5Dopen(group_id, name.c_str(), H5P_DEFAULT);
    
    if(id<0) 
      throw std::logic_error(__FUNCTION__);
  }
    
  hdf5_dataset::~hdf5_dataset()
  {
    if(H5Dclose(id)<0) 
      throw std::logic_error(__FUNCTION__);
  }

  hid_t hdf5_dataset::get_id()
  {
    return id;
  }

}

#endif
