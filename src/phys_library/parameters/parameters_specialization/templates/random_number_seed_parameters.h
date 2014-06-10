//-*-C++-*-

#ifndef RANDOM_NUMBER_SEED_PARAMETERS_H
#define RANDOM_NUMBER_SEED_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class random_number_seed_parameters 
{
public:

  random_number_seed_parameters();
  ~random_number_seed_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size(const concurrency_type& concurrency) const;

  template<class concurrency_type>
  void pack           (const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         (const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

  template<class stream_type>
  void to_JSON(stream_type& ss, bool is_end=false);
  
  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

/******************************************
 ***        DATA                        ***
 ******************************************/

//   int get_seed();

private:

  int seed;

};

random_number_seed_parameters::random_number_seed_parameters():
  seed(985456376)
{}

random_number_seed_parameters::~random_number_seed_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int random_number_seed_parameters::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(seed);

  return buffer_size;
}

template<class concurrency_type>
void random_number_seed_parameters::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, seed);
}

template<class concurrency_type>
void random_number_seed_parameters::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, seed);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void random_number_seed_parameters::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"random-number-seed\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "seed", seed, true);
  
  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";
}
  
template<class JSON_reader_type>
void random_number_seed_parameters::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  
  try
    {
      const JsonAccessor control(reader["random-number-seed"]);
      seed <= control["seed"];
    }
  catch(const std::exception& r_e)
    {
      cout << "\n\n\trandom-number-seed --> " << seed << "\n\n";
    }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

// int random_number_seed_parameters::get_seed()
// {
//   return seed;
// }

#endif
