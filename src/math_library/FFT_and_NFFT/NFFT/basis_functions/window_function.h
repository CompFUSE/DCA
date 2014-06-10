//-*-C++-*-

#ifndef MATH_ALGORITHMS_NFFT_WINDOW_FUNCTION_H
#define MATH_ALGORITHMS_NFFT_WINDOW_FUNCTION_H

namespace MATH_ALGORITHMS
{
  namespace NFFT
  {
    enum NFFT_WINDOW_FUNCTION_NAMES {GAUSSIAN, KAISER_BESSEL};

    template<typename dnfft_type>
    struct nfft_window_function
    {
      typedef typename dnfft_type::scalar_type scalar_type;

    public:

      nfft_window_function();

      ~nfft_window_function();

      void initialize(NFFT_WINDOW_FUNCTION_NAMES NAME,
		      int                        n,
		      int                        m,
		      scalar_type                sigma);

      scalar_type phi_t  (scalar_type t);
      scalar_type d_phi_t(scalar_type t);

      scalar_type phi_wn (int wn);

    private:

      NFFT_WINDOW_FUNCTION_NAMES NAME;

      int    n;
      int    m;

      scalar_type sigma;
    };

    template<typename dnfft_type>
    nfft_window_function<dnfft_type>::nfft_window_function():
    {}

    template<typename dnfft_type>
    nfft_window_function<dnfft_type>::~nfft_window_function()
    {}

    template<typename dnfft_type>
    void nfft_window_function<dnfft_type>::set_type(NFFT_WINDOW_FUNCTION_NAMES NEW_NAME,
						    int                        max_frequencies,
						    int                        over_sampling,
						    scalar_type                sigma_value)
    {
      NAME = NEW_NAME;

      switch(NAME)
	{
	case KAISER_BESSEL:
	  kaiser_bessel_function::n = max_frequencies;
	  kaiser_bessel_function::m = over_sampling;
	  
	  kaiser_bessel_function::sigma = sigma_values;
	  break;
	  
	default:
	  throw std::logic_error(__FUNCTION__);
	}
    }

    template<typename dnfft_type>
    void nfft_window_function<dnfft_type>::phi_t(scalar_type t)
    {
      switch(NAME)
	{
	case KAISER_BESSEL:
	  return kaiser_bessel_function::phi_t(t);
	  break;

	default:
	  throw std::logic_error(__FUNCTION__);
	}
    }

    template<typename dnfft_type>
    void nfft_window_function<dnfft_type>::d_phi_t(scalar_type t)
    {
      switch(NAME)
	{
	case KAISER_BESSEL:
	  return kaiser_bessel_function::d_phi_t(t);
	  break;

	default:
	  throw std::logic_error(__FUNCTION__);
	}
    }

    template<typename dnfft_type>
    void nfft_window_function<dnfft_type>::phi_wn(int wn)
    {
      switch(NAME)
	{
	case KAISER_BESSEL:
	  return kaiser_bessel_function::phi_wn(t);
	  break;

	default:
	  throw std::logic_error(__FUNCTION__);
	}

    }

  }

}

#endif
