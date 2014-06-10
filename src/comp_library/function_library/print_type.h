//-*-C++-*-

#ifndef PRINT_SCALAR_TYPE_H
#define PRINT_SCALAR_TYPE_H

/*
template<typename scalartype>
class print_scalar_type
{
public:

  static void execute(std::ofstream& ss, scalartype& f)
  {
    ss << real(do_cast<std::complex<double> >::execute(f));
    ss << ",";
    ss << imag(do_cast<std::complex<double> >::execute(f));
  }

};

template<typename scalartype>
class print_scalar_type<std::vector<scalartype> >
{
public:

  static void execute(std::ofstream& ss, std::vector<scalartype>& f)
  {
    for(size_t l=0; l<f.size(); l++){
      ss << do_cast<double>::execute(f[l]);

      if(l<f.size()-1)
	ss << ",";
    }
  }

};
*/

#endif
