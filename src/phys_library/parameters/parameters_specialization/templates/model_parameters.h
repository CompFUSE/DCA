//-*-C++-*-

#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
template <class model_t>
class model_parameters {
public:
  model_parameters();
  ~model_parameters();

  /******************************************
   ***        CONCURRENCY                 ***
   ******************************************/

  template <class concurrency_type>
  int get_buffer_size(const concurrency_type& concurrency) const;

  template <class concurrency_type>
  void pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template <class concurrency_type>
  void unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  /******************************************
   ***        READ/WRITE                  ***
   ******************************************/

  template <class stream_type>
  void to_JSON(stream_type& ss, bool is_end = false);

  template <class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template <class read_write_type>
  void read_write(read_write_type& read_write_obj);

  /******************************************
   ***        DATA                        ***
   ******************************************/
  double get_t();
  double get_t_prime();
  double get_V();
  double get_V_prime();
  double get_U();
};

#endif
