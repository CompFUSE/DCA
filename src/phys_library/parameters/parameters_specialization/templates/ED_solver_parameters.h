//-*-C++-*-

#ifndef ED_SOLVER_PARAMETERS_H
#define ED_SOLVER_PARAMETERS_H

/*!
 *   \author Peter Staar
 */
class ED_solver_parameters
{
public:

  ED_solver_parameters();
  ~ED_solver_parameters();

  /******************************************
   ***        CONCURRENCY                 ***
   ******************************************/

  template<class concurrency_type>
  int  get_buffer_size( concurrency_type& concurrency);

  template<class concurrency_type>
  void pack           ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  /******************************************
   ***        READ/WRITE                  ***
   ******************************************/

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

  /******************************************
   ***        DATA                        ***
   ******************************************/

  double      get_eigenvalue_cut_off();

  bool        do_orthogonality_check();

  std::string get_symmetries();

  std::string get_ED_method();

  bool        do_sector_check();

  int         get_occupation();

  int         get_magnetization();

  int         get_n_0();
  int         get_Sz_0();
  int         get_n_1();
  int         get_Sz_1();
  int         get_nu();

private:

  double      eigenvalue_cut_off;

  std::string check_orthogonality_of_states;

  std::string symmetries;

  std::string ED_method;

  std::string check_sector;

  int         occupation;

  int         magnetization;

  int         n_0;
  int         Sz_0;
  int         n_1;
  int         Sz_1;
  int         nu;

};

ED_solver_parameters::ED_solver_parameters():
  eigenvalue_cut_off(1.e-6),

  check_orthogonality_of_states("false"),

  symmetries("default"),

  ED_method("default"),

  check_sector("false"),

  occupation(0),

  magnetization(0)

{}

ED_solver_parameters::~ED_solver_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int ED_solver_parameters::get_buffer_size(concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(eigenvalue_cut_off);
  buffer_size += concurrency.get_buffer_size(ED_method);

  return buffer_size;
}

template<class concurrency_type>
void ED_solver_parameters::pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, eigenvalue_cut_off);
  concurrency.pack(buffer, buffer_size, position, ED_method);
}

template<class concurrency_type>
void ED_solver_parameters::unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, eigenvalue_cut_off);
  concurrency.unpack(buffer, buffer_size, position, ED_method);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class read_write_type>
void ED_solver_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("ED-solver-parameters");

      try { read_write_obj.execute("eigenvalue-cut-off", eigenvalue_cut_off); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("check_orthogonality_of_states", check_orthogonality_of_states); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("symmetries", symmetries); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("ED_method", ED_method); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("check_sector", check_sector); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("occupation", occupation); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("magnetization", magnetization); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("n_0", n_0); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("Sz_0", Sz_0); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("n_1", n_1); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("Sz_1", Sz_1); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("nu", nu); } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e)
    {}
}

double ED_solver_parameters::get_eigenvalue_cut_off()
{
  return eigenvalue_cut_off;
}

bool ED_solver_parameters::do_orthogonality_check()
{
  if (check_orthogonality_of_states == "true")
    return true;
  else
    return false;
}

std::string ED_solver_parameters::get_symmetries()
{
  return symmetries;
}

std::string ED_solver_parameters::get_ED_method()
{
  return ED_method;
}

bool ED_solver_parameters::do_sector_check()
{
  if (check_sector == "true")
    return true;
  else
    return false;
}

int ED_solver_parameters::get_occupation()
{
  return occupation;
}

int ED_solver_parameters::get_magnetization()
{
  return magnetization;
}

int ED_solver_parameters::get_n_0()
{
  return n_0;
}
int ED_solver_parameters::get_Sz_0()
{
  return Sz_0;
}
int ED_solver_parameters::get_n_1()
{
  return n_1;
}
int ED_solver_parameters::get_Sz_1()
{
  return Sz_1;
}
int ED_solver_parameters::get_nu()
{
  return nu;
}

#endif
