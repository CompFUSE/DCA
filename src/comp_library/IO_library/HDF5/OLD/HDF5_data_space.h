//-*-C++-*-

#ifndef HDF5_DATASPACE_HEADER_H
#define HDF5_DATASPACE_HEADER_H

namespace IO
{
  class hdf5_dataspace
  {
  public:
 
    hdf5_dataspace(const std::vector<int> dims);
    ~hdf5_dataspace();

    hid_t get_id();

  private:
    
    hid_t id;
  };

  hdf5_dataspace::hdf5_dataspace(const std::vector<int> dims)
  {
    std::vector<hsize_t> current_dims(dims.size());
    
    for(int i=0; i<(int)dims.size(); i++) 
      current_dims[dims.size() - i - 1] = dims[i];
    
    id = H5Screate_simple((int)dims.size(), &current_dims[0], NULL);
    
    if(id<0) 
      throw std::logic_error(__FUNCTION__);
  }

  hdf5_dataspace::~hdf5_dataspace()
  {
    if(H5Sclose(id) < 0) 
      throw std::logic_error(__FUNCTION__);
  }

  hid_t hdf5_dataspace::get_id()
  {
    return id;
  }

}

#endif
