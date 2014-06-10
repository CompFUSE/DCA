//-*-C++-*-

#ifndef NUMERICAL_ERROR_DOMAIN_H
#define NUMERICAL_ERROR_DOMAIN_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \author Peter Staar
     */
    class numerical_error_domain
    {
    public:

      typedef double                 element_type;
      typedef numerical_error_domain this_type;

    public:

      static int                       get_size();
      static std::vector<element_type>& get_elements();

      template<class stream_type>
      static void to_JSON(stream_type& ss);

    private:

      static std::vector<element_type>& initialize_elements();
    };

    int numerical_error_domain::get_size()
    {
      return int(get_elements().size());
    }

    std::vector<double>& numerical_error_domain::get_elements()
    {
      static std::vector<element_type>& v = initialize_elements();
      return v;
    }

    std::vector<double>& numerical_error_domain::initialize_elements()
    {
      static std::vector<element_type> v(0);

      for(int i=-16; i<0; i++)
        for(double j=1; j<10; j+=1.)
          v.push_back(j*std::pow(10., i));

      return v;
    }

    template<class stream_type>
    void numerical_error_domain::to_JSON(stream_type& ss)
    {
      ss << "\"numerical_error_domain\" : [\n";

      for(int i=0; i<get_size(); i++)
        if(i == get_size()-1)
          ss << get_elements()[i] << "\n";
        else
          ss << get_elements()[i] << ",\n";

      ss << "]\n";
    }

  }

}

#endif
