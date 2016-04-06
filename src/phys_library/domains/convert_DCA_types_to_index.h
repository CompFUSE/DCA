//-*-C++-*-
  /*!
   *  \author: Peter Staar
   */
#ifndef CONVERT_DCA_TYPES_TO_INDEX_H
#define CONVERT_DCA_TYPES_TO_INDEX_H

namespace QMC
{

  template<typename target, class whatever>
  class convert{};
  template<>
  class convert<int, dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > >
  {
    typedef dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > nu;

  public:

    static int spin_orbital(int band, e_spin_states_type e_spin);

  private:
    static FUNC_LIB::function<int, nu>& intitialize_spin_orbital();
  };

  int convert<int, dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > >::spin_orbital(int band,
                                                                                                   e_spin_states_type e_spin)
  {
    typedef dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > nu;
    typedef nu::this_type parameter_typelist;

    static FUNC_LIB::function<int, nu>& spo_function = intitialize_spin_orbital();

    if(std::is_same<TypeAt<0,parameter_typelist>::type, electron_band_domain>::value)
      {
        return spo_function(band, electron_spin_domain::to_coordinate(e_spin));
        //return spo_function(band, do_cast<int>::execute(e_spin));
      }
    else
      {
        return spo_function(electron_spin_domain::to_coordinate(e_spin), band);
        //return spo_function(do_cast<int>::execute(e_spin), band);
      }
  }

  FUNC_LIB::function<int, dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > >&
  convert<int, dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > >::intitialize_spin_orbital()
  {
    typedef dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > nu;
    static FUNC_LIB::function<int, nu> spo_function;

    for(int i=0; i<spo_function.size(); i++)
      spo_function(i) = i;

    return spo_function;
  }

  /*!
   *  \author: Peter Staar
   */
  template<>
  class convert<dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> >, int>
  {
    typedef dmn_0<electron_band_domain> b;
    typedef dmn_0<electron_spin_domain> s;

    typedef dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> > nu;

  public:

    static std::pair<int,int> from_spin_orbital(int spo);
  };

  std::pair<int,int> convert<dmn_2<dmn_0<electron_band_domain>,dmn_0<electron_spin_domain> >, int>::from_spin_orbital(int spo)
  {
    std::pair<int,int> tmp;
    tmp.first  = spo % b::dmn_size();
    tmp.second = (spo-tmp.first) % s::dmn_size();
    return tmp;
  }

}

#endif
