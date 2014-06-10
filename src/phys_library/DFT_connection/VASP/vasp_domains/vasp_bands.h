//-*-C++-*-

#ifndef VASP_BAND_DOMAIN_H
#define VASP_BAND_DOMAIN_H

namespace DFT
{
  namespace VASP
  {
    /*!
     *  \author Peter Staar
     */
    class vasp_band_domain
    {
    public:

      typedef int element_type;

    public:

      static int&                       get_size();
      static std::string                get_name();

      static std::vector<element_type>& get_elements();

      template<IO::FORMAT DATA_FORMAT>
      static void read(IO::reader<DATA_FORMAT>& reader);

      template<IO::FORMAT DATA_FORMAT>
      static void write(IO::writer<DATA_FORMAT>& writer);

      template<typename parameters_type>
      static void initialize(parameters_type& parameters);

      template<class stream_type>
      static void to_JSON(stream_type& ss);

    private:

      static std::vector<element_type>& initialize_elements();
    };

    int& vasp_band_domain::get_size()
    {
      static int size = 0;
      return size;
    }

    std::string vasp_band_domain::get_name()
    {
      static std::string name = "vasp-band-domain";
      return name;
    }

    std::vector<vasp_band_domain::element_type>& vasp_band_domain::get_elements()
    {
      static std::vector<element_type> elements(get_size());
      return elements;
    }

    template<IO::FORMAT DATA_FORMAT>
    void vasp_band_domain::write(IO::writer<DATA_FORMAT>& writer)
    {
      writer.open_group(get_name());
      writer.execute(get_elements());
      writer.close_group();
    }

    template<typename parameters_type>
    void vasp_band_domain::initialize(parameters_type& parameters)
    {
      get_size() = parameters.get_nb_vasp_bands();

      for(size_t i=0; i<get_elements().size(); ++i){
        get_elements()[i] = i;
      }
    }

  }

}

#endif
