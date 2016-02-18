//-*-C++-*-

#ifndef PARALLELIZATION_LIBRARY_MPI_H
#define PARALLELIZATION_LIBRARY_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class parallelization<MPI_LIBRARY> : public print_on_shell_interface<MPI_LIBRARY>,
                                       public packing_interface       <MPI_LIBRARY>,
                                       public RNG_interface           <MPI_LIBRARY>,
                                       public collective_min_interface<MPI_LIBRARY>,
                                       public collective_max_interface<MPI_LIBRARY>,
                                       public collective_sum_interface<MPI_LIBRARY>
  {
  public:

    parallelization(int argc, char *argv[]);
    ~parallelization();

    int id();

    int number_of_processors();

    int first();
    int last ();

    void set_seed(long seed);

    long get_seed();

    template<typename object_type>
    bool broadcast(object_type& object, int root_id=0);
    
    template<typename object_type>
    bool broadcast_object(object_type& object, int root_id=0);

    template<typename domain_type>
    std::pair<int, int> get_bounds(domain_type& dmn);

  private:

    processor_grouping<MPI_LIBRARY> group;

    long SEED;
  };

  parallelization<MPI_LIBRARY>::parallelization(int argc, char *argv[]):
    print_on_shell_interface<MPI_LIBRARY>(group),
    packing_interface       <MPI_LIBRARY>(group),
    RNG_interface           <MPI_LIBRARY>(group),
    collective_min_interface<MPI_LIBRARY>(group),
    collective_max_interface<MPI_LIBRARY>(group),
    collective_sum_interface<MPI_LIBRARY>(group),

    SEED(0)
  {
    MPI_Init(&argc, &argv);

    group.set();

    //   {
    //     std::stringstream ss;
    //     ss << get_id() << "\t" << get_Nr_threads() << endl;
    //     std::cout << ss.str();
    //   }

    RNG_interface<MPI_LIBRARY>::initialize();
  }

  parallelization<MPI_LIBRARY>::~parallelization()
  {
    MPI_Finalize();
  }

  int parallelization<MPI_LIBRARY>::first()
  {
    return group.first();
  }
  
  int parallelization<MPI_LIBRARY>::last()
  {
    return group.last();
  } 

  int parallelization<MPI_LIBRARY>::id()
  {
    return group.get_id();
  }

  int parallelization<MPI_LIBRARY>::number_of_processors()
  {
    assert(group.get_Nr_threads()>-1);
    return group.get_Nr_threads();
  }

  long parallelization<MPI_LIBRARY>::get_seed()
  {
    return SEED;
  }

  void parallelization<MPI_LIBRARY>::set_seed(long seed)
  {
    srand(seed);
    
    for(int l=0; l<=id(); l++)
      SEED = rand();
  }


  template<typename object_type>
  bool parallelization<MPI_LIBRARY>::broadcast(object_type& object, int root_id)
  {
    assert(root_id>-1 and root_id<number_of_processors());

    int position = 0;

    if(id() == root_id)
      {
        int bufferSize = this->get_buffer_size(object);

        MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, group.get());

        int* buffer = new int[bufferSize];

        this->pack(buffer, bufferSize, position, object);

        MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, group.get());

        delete [] buffer;
      }
    else
      {
        int bufferSize(0);

        MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, group.get());

        int* buffer = new int[bufferSize];

        // receive packed message
        MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, group.get());

        this->unpack(buffer, bufferSize, position, object);

        delete [] buffer;
      }

    return true;
  }

  template<typename object_type>
  bool parallelization<MPI_LIBRARY>::broadcast_object(object_type& object, int root_id)
  {
    assert(root_id>-1 and root_id<number_of_processors());

    int buffer_size = 0;

    if(id() == root_id)
      {
        buffer_size = object.get_buffer_size(*this);

        MPI_Bcast(&buffer_size, 1          , MPI_INT   , root_id, group.get());

        int* buffer = new int[buffer_size];

        int off_set = 0;
	object.pack(*this, buffer, buffer_size, off_set);
        //object.pack_or_unpack(true, *this, buffer, buffer_size, off_set);

        MPI_Bcast(buffer     , buffer_size, MPI_PACKED, root_id, group.get());

        delete [] buffer;
      }
    else
      {
        MPI_Bcast(&buffer_size, 1          , MPI_INT   , root_id, group.get());

        int* buffer = new int[buffer_size];

        MPI_Bcast(buffer      , buffer_size, MPI_PACKED, root_id, group.get());

        int off_set = 0;
	object.unpack(*this, buffer, buffer_size, off_set);
        //object.pack_or_unpack(false, *this, buffer, buffer_size, off_set);

        delete [] buffer;
      }

    return true;
  }

  template<typename domain_type>
  std::pair<int, int> parallelization<MPI_LIBRARY>::get_bounds(domain_type& dmn)
  {
    long long size = static_cast<long long>(dmn.get_size());

    long long bounds_first, bounds_second;

    long long CPU_id = static_cast<long long>(id());
    long long np     = static_cast<long long>(number_of_processors());

    if(np < size)
      {
        bounds_first  = ( CPU_id   *size)/np;
        bounds_second = ((CPU_id+1)*size)/np;
      }
    else
      {
        if(CPU_id < size){
          bounds_first  = CPU_id;
          bounds_second = CPU_id+1;
        }
        else{
          bounds_first  = -1;
          bounds_second = -1;
        }
      }

    std::pair<int, int> bounds(static_cast<int>(bounds_first), static_cast<int>(bounds_second));

    if(!( (bounds.first==-1 && bounds.second==-1) || (bounds.first>=0 && bounds.second<=dmn.get_size() && bounds.first<bounds.second) ) )
      {
        std::cout << "error in " << __PRETTY_FUNCTION__ << "\n\n";
        std::cout << "CPU-id :: " << CPU_id << "\n";
        std::cout << "np     :: " << np     << "\n";

        std::cout << "bounds.first  :: " << bounds.first  << "\n";
        std::cout << "bounds.second :: " << bounds.second << "\n";

        throw std::logic_error(__FUNCTION__);
      }

    return bounds;
  }

}

#endif
