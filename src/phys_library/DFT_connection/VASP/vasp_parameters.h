//-*-C++-*-

#ifndef DFT_VASP_PARAMETERS_H
#define DFT_VASP_PARAMETERS_H

namespace DFT
{
  namespace VASP
  {

    /*!
     *   \ingroup  PARAMETERS
     *
     *   \author   Peter Staar, Long Zhang
     *   \brief    ...
     */
    template<class concurrency_t>
    class parameters
    {
#include "type_definitions.h"

      typedef concurrency_t concurrency_type;

    public:

      parameters(std::string       version_stamp,
                 concurrency_type& concurrency_obj);

      ~parameters();

      template<IO::FORMAT DATA_FORMAT>
      void read(IO::reader<DATA_FORMAT>& reader);

      template<IO::FORMAT DATA_FORMAT>
      void write(IO::writer<DATA_FORMAT>& writer);

      void read_input(string file_name);

      void write_input(string file_name);

      template<typename read_write_t>
      void read_write(read_write_t& read_write_obj);

      void update_domains();

      int  get_buffer_size(concurrency_type& concurrency);
      void pack           (concurrency_type& concurrency, int* buffer, int buffer_size, int& position);
      void unpack         (concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

      // get_functions
      concurrency_type& get_concurrency();

      double get_epsilon() { return 1.e-6; }

      std::string get_PROCAR_file()         { return PROCAR_file_name;}
      std::string get_t_ij_file()           { return t_ij_file_name;}
      std::string get_band_structure_file() { return band_structure_file_name;}

      int get_nb_dmft_bands()    { return nCorrBands;}
      int get_nb_dmft_orbitals() { return CorrOrbiIndex.size()/*nCorrBands*/;}

      int get_nb_vasp_bands()    { return nBands;}
      int get_nb_vasp_orbitals() { return nOrbitals;}

      int get_iid() { return iid; }

      int get_nBands  () { return nBands; }
      int get_nKpoints() { return k_vasp::dmn_size(); }

      int get_vasp_orbitals() { return nOrbitals; }

      double get_Fermi_level() { return EFermiDFT;}

      int get_lower_index_of_chosen_bands() { return corrband_min; }
      int get_upper_index_of_chosen_bands() { return corrband_max; }

      std::vector<int> get_correlated_orbitals() { return CorrOrbiIndex; }

      int get_nr() { return nr; }

      std::string&                       get_coordinate_type() { return coordinate_type; }
      std::vector<std::string>&          get_coordinate_names() { return coordinate_names; }
      std::vector<std::vector<double> >& get_Brillouin_zone_vectors() { return coordinate_vectors; }


    private:

      //       template<typename read_write_type>
      //       void read_write(read_write_type& obj);

    private:

      concurrency_type& concurrency;

      std::string PROCAR_file_name;
      std::string t_ij_file_name;
      std::string band_structure_file_name;

      double Ewindow_up;      // for finding the desired bands
      double Ewindow_dn;      // for finding the desired bands

      int    nCorrBands;         // number of chosen bands, need this for checking # of lines
      int  corrband_min;        // lower index of chosen bands ( range [0,nBands-1] )
      int  corrband_max;        // upper index of chosen bands

      int      nKpoints;       // look at the top of PROCAR, or KPOINTS
      int        nBands;        // look at the top of PROCAR
      int         nIons;         // look at the top of PROCAR
      double      alatt;     // look at POSCAR
      double  EFermiDFT;    // grep "E-fermi" OUTCAR

      int            nr;         // range of neighbors cells, for calculating tij

      int           iid;         // choice of correlation site(atom), 2 for SrVO3, 1 for NiO

      int nOrbitals;

      std::vector<int> K_grid;      // need to m
      std::vector<int> CorrOrbiIndex;      // need to manually assign value to this

      double* a_lattice;

      std::string                       coordinate_type;

      std::vector<std::string>          coordinate_names;
      std::vector<std::vector<double> > coordinate_vectors;
    };

    template<class concurrency_t>
    parameters<concurrency_t>::parameters(std::string       version_stamp,
                                          concurrency_type& concurrency_ref):
      concurrency(concurrency_ref),

      PROCAR_file_name        ("default-PROCAR-file-name"),
      t_ij_file_name          ("t-ij.txt"),
      band_structure_file_name("band-structure.ps"),

      Ewindow_up(-1),      // for finding the desired bands
      Ewindow_dn( 1),      // for finding the desired bands

      nCorrBands(1),         // number of chosen bands, need this for checking # of lines
      corrband_min(0),        // lower index of chosen bands ( range [0,nBands-1] )
      corrband_max(1),        // upper index of chosen bands

      nKpoints(64),       // look at the top of PROCAR, or KPOINTS
      nBands(1),        // look at the top of PROCAR
      nIons(1),         // look at the top of PROCAR
      alatt(1.),     // look at POSCAR
      EFermiDFT(0.),    // grep "E-fermi" OUTCAR

      nr(5),         // range of neighbors cells, for calculating tij

      iid(2),         // choice of correlation site(atom), 2 for SrVO3, 1 for NiO

      nOrbitals(1),

      K_grid(3, 4),      // need to m
      CorrOrbiIndex(1, 0),      // need to manually assign value to this

      a_lattice(NULL),

      coordinate_type("default")
    {
      a_lattice = new double[3*3];

      for(int j=0; j<3; j++)
        for(int i=0; i<3; i++)
          a_lattice[i+3*j] = i==j? 1 : 0;
    }

    template<class concurrency_t>
    parameters<concurrency_t>::~parameters()
    {
      delete [] a_lattice;
    }

    template<class concurrency_t>
    void parameters<concurrency_t>::read_input(string file_name)
    {
      //if(concurrency_obj.id() == concurrency_obj.first())
      {
        //file_names_parameters::get_input_file_name() = filename;

        {
          IO::reader<IO::JSON> read_obj;

          read_obj.open_file(file_name);

          read_write(read_obj);

          read_obj.close_file();
        }

      }

      //       {
      //        for(int d=0; d<3; d++)
      //          cout << d << "\t" << K_grid[d] << "\n";
      //        cout << "\n";
      //       }
    }

    template<class concurrency_t>
    void parameters<concurrency_t>::write_input(string file_name)
    {
      IO::writer<IO::JSON> write_obj;

      write_obj.open_file(file_name);

      read_write(write_obj);

      write_obj.close_file();
    }

    template<class concurrency_t>
    template<typename read_write_type>
    void parameters<concurrency_t>::read_write(read_write_type& obj)
    {
      {
        obj.open_group("file-parameters");

        try { obj.execute("t_ij-file"          , t_ij_file_name); } catch(const std::exception& r_e) {}
        try { obj.execute("band-structure-file", band_structure_file_name); } catch(const std::exception& r_e) {}

        obj.close_group();
      }

      {
        obj.open_group("DMFT-parameters");

        try { obj.execute("cut-off-radius", nr); } catch(const std::exception& r_e) {}

        try { obj.execute("iid", iid); } catch(const std::exception& r_e) {}

        try { obj.execute("w-min", Ewindow_up); } catch(const std::exception& r_e) {}
        try { obj.execute("w-max", Ewindow_dn); } catch(const std::exception& r_e) {}

        try { obj.execute("orbitals", CorrOrbiIndex); } catch(const std::exception& r_e) {}

        try { obj.execute("lower-index-of-chosen-bands", corrband_min); } catch(const std::exception& r_e) {}
        try { obj.execute("upper-index-of-chosen-bands", corrband_max); } catch(const std::exception& r_e) {}

        obj.close_group();
      }

      {
        obj.open_group("VASP-parameters");

        try { obj.execute("PROCAR-file", PROCAR_file_name); } catch(const std::exception& r_e) {}
        try { obj.execute("K-grid"     , K_grid);         } catch(const std::exception& r_e) {}

        try { obj.execute("# bands"    , nBands);           } catch(const std::exception& r_e) {}
        try { obj.execute("# orbitals" , nOrbitals);           } catch(const std::exception& r_e) {}

        try { obj.execute("# ions"     , nIons);            } catch(const std::exception& r_e) {}

        //try { obj.execute("# a-lattice"  , alatt);          } catch(const std::exception& r_e) {}
        try { obj.execute("Fermi-level", EFermiDFT);      } catch(const std::exception& r_e) {}

        obj.close_group();
      }

      try
        {
          obj.open_group("band-structure-cut");

          try { obj.execute("coordinate-type"       , coordinate_type);    } catch(const std::exception& r_e) {}
          try { obj.execute("Brillouin-zone-names"  , coordinate_names);   } catch(const std::exception& r_e) {}
          try { obj.execute("Brillouin-zone-vectors", coordinate_vectors); } catch(const std::exception& r_e) {}

          obj.close_group();
        }
      catch(const std::exception& r_e)
        {}

      {
        if(corrband_max<corrband_min)
          {
            cout << "\n\n upper-index-of-chosen-bands < lower-index-of-chosen-bands is not allowed!!!";
            throw std::logic_error(__FUNCTION__);
          }

        nCorrBands = corrband_max-corrband_min+1;

        if(nCorrBands<CorrOrbiIndex.size())
          {
            cout << "\n\n (upper-index-of-chosen-bands - lower-index-of-chosen-bands + 1) < size(orbitals) is not allowed!!!";
            throw std::logic_error(__FUNCTION__);
          }

      }
    }

    template<class concurrency_t>
    void  parameters<concurrency_t>::update_domains()
    {
      b_dmft::parameter_type::initialize(*this);
      o_dmft::parameter_type::initialize(*this);

      b_vasp::parameter_type::initialize(*this);
      o_vasp::parameter_type::initialize(*this);

      cluster_domain_initializer<r_vasp>::execute(a_lattice, K_grid);
    }

  }

}

#endif

