//-*-C++-*-

#ifndef SET_TO_ZERO
#define SET_TO_ZERO

/*! \class  set_to_zero
 *  \author Peter Staar, Bart Ydens
 *  \brief  This class will assure that all functions with a scalartype of float, double, std::complex<float> and std::complex<double> are initialized to zero. Functions with a different type will not be initialized!
 *  \date 2011
 *  \version 1.0
 */

struct set_to_zero
{
public:

  template<class whatever_t>
  static void execute(whatever_t& whatever);
};

template<class whatever_t>
void set_to_zero::execute(whatever_t& whatever)
{}

template<>
void set_to_zero::execute(short& whatever)
{
  whatever = 0;
}

template<>
void set_to_zero::execute(int& whatever)
{
  whatever = 0;
}

template<>
void set_to_zero::execute(long& whatever)
{
  whatever = 0;
}

template<>
void set_to_zero::execute(size_t& whatever)
{
  whatever = 0;
}

template<>
void set_to_zero::execute(float& whatever)
{
  whatever = 0.;
}

template<>
void set_to_zero::execute(double& whatever)
{
  whatever = 0.;
}

template<>
void set_to_zero::execute(std::complex<float>& whatever)
{
  whatever = 0.;
}

template<>
void set_to_zero::execute(std::complex<double>& whatever)
{
  whatever = 0.;
}

#endif
