// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Long Zhang
//
// VASP reader.

#ifndef DCA_PHYS_DFT_CONNECTION_VASP_READER_HPP
#define DCA_PHYS_DFT_CONNECTION_VASP_READER_HPP

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dft_connection/vasp/data.hpp"
#include "dca/phys/dft_connection/vasp/vasp_domains/dmft_band_domain.hpp"
#include "dca/phys/dft_connection/vasp/vasp_domains/dmft_orbital_domain.hpp"
#include "dca/phys/dft_connection/vasp/vasp_domains/vasp_band_domain.hpp"
#include "dca/phys/dft_connection/vasp/vasp_domains/vasp_orbital_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"

namespace dca {
namespace phys {
namespace dft {
namespace vasp {
// dca::phys::dft::vasp::

template <class parameters_type>
class reader {
public:
  using k_vasp =
      func::dmn_0<domains::cluster_domain<double, 3, domains::VASP_LATTICE, domains::MOMENTUM_SPACE,
                                          domains::PARALLELLEPIPEDUM>>;

  using b_dmft = func::dmn_0<dmft_band_domain>;
  using o_dmft = func::dmn_0<dmft_orbital_domain>;

  using b_vasp = func::dmn_0<vasp_band_domain>;
  using o_vasp = func::dmn_0<vasp_orbital_domain>;

public:
  reader(parameters_type& parameters, data<parameters_type>& data);

  ~reader();

  void execute();

private:
  parameters_type& parameters;

  data<parameters_type>& vasp_data;

  func::function<double, k_vasp>& kx;
  func::function<double, k_vasp>& ky;
  func::function<double, k_vasp>& kz;

  func::function<double, func::dmn_variadic<k_vasp, b_vasp>>& band_e;
  func::function<double, func::dmn_variadic<k_vasp, b_vasp>>& band_o;

  func::function<double, func::dmn_variadic<k_vasp, b_vasp, o_vasp>>& proj_magni;  // proj_magni[nKpoints][nBands][nOrbi]
  func::function<double, func::dmn_variadic<k_vasp, b_vasp, o_vasp>>& proj_phaRe;  // proj_phaRe[nKpoints][nBands][nOrbi]
  func::function<double, func::dmn_variadic<k_vasp, b_vasp, o_vasp>>& proj_phaIm;  // proj_phaIm[nKpoints][nBands][nOrbi]

  func::function<std::complex<double>, func::dmn_variadic<k_vasp, b_dmft, o_dmft>>& projector;

  //       func::function<double, func::dmn_variadic<k_vasp, b_dmft, o_dmft> >& projection_re;
  //       //
  //       projection_Re[nKpoints][nCorrBands][nCorrOrbis]
  //       func::function<double, func::dmn_variadic<k_vasp, b_dmft, o_dmft> >& projection_im;
  //       //
  //       projection_Im[nKpoints][nCorrBands][nCorrOrbis]

  //       func::function<double, func::dmn_variadic<k_vasp, b_dmft, b_dmft> >& BlochHami; //
  //       BlochHami[nKpoints][nCorrBands][nCorrBands]
  //       func::function<double, func::dmn_variadic<k_vasp, b_dmft, b_dmft> >& BlochHami_re;//
  //       dfBlochHamiRe
  //       [nKpoints][nCorrOrbis][nCorrOrbis]
  //       func::function<double, func::dmn_variadic<k_vasp, b_dmft, b_dmft> >& BlochHami_im;//
  //       dfBlochHamiRe
  //       [nKpoints][nCorrOrbis][nCorrOrbis]

  /*!
   *   extra data for the parsing ...
   */

  std::string target1;
  std::string target2;
  std::string target3;
  std::string target4;
  std::string target5;

  std::string tagk1;
  std::string tagk2;
  std::string tagb1;
  std::string tagb2;

  double s_R, py_R, pz_R, px_R, dxy_R, dyz_R, dz2_R, dxz_R, dx2_R;
  double s_I, py_I, pz_I, px_I, dxy_I, dyz_I, dz2_I, dxz_I, dx2_I;

  double xx, yy;

  int counter1, counter2, counter3, counter4;
  std::string oneline, dummy, dummy2;

  int pos0, pos1, pos2, pos3;
  int ndp1, ndp2, ndp3;
  int nstar1, nstar2, nstar3;
  std::vector<int> head;
  std::vector<int> tail;
  std::vector<int> pdp;

  // const double PI = 3.14159265358979323846;

  std::vector<std::string> useful_lines;
};

template <class parameters_type>
reader<parameters_type>::reader(parameters_type& parameters_ref, data<parameters_type>& vasp_data_ref)
    : parameters(parameters_ref),
      vasp_data(vasp_data_ref),

      kx(vasp_data.kx),
      ky(vasp_data.ky),
      kz(vasp_data.kz),

      band_e(vasp_data.band_e),
      band_o(vasp_data.band_o),

      proj_magni(vasp_data.proj_magni),  // proj_magni[nKpoints][nBands][nOrbi]
      proj_phaRe(vasp_data.proj_phaRe),  // proj_phaRe[nKpoints][nBands][nOrbi]
      proj_phaIm(vasp_data.proj_phaIm),  // proj_phaIm[nKpoints][nBands][nOrbi]

      projector(vasp_data.projector),

      //       projection_re(vasp_data.projection_re), //
      //       projection_Re[nKpoints][nCorrBands][nCorrOrbis]
      //       projection_im(vasp_data.projection_im), //
      //       projection_Im[nKpoints][nCorrBands][nCorrOrbis]

      //       BlochHami(vasp_data.BlochHami), // BlochHami[nKpoints][nCorrBands][nCorrBands]

      //       BlochHami_re(vasp_data.BlochHami_re),// dfBlochHamiRe
      //       [nKpoints][nCorrOrbis][nCorrOrbis]
      //       BlochHami_im(vasp_data.BlochHami_im),// dfBlochHamiRe
      //       [nKpoints][nCorrOrbis][nCorrOrbis]

      target1("ion"),  // tag the lines that are not useful
      target2("tot"),  // only my choices
      target3("PROCAR lm decomposed"),
      target4("# of bands:"),
      target5("# of k-points:"),

      tagk1("k-point"),  // tag k-point line
      tagk2("weight"),   // tag k-point line
      tagb1("band"),     // tag band energy line
      tagb2("occ.")      // tag band energy line
{
  std::cout << "\n PART I starts. \n";
}

template <class parameters_type>
reader<parameters_type>::~reader() {
  std::cout << "\n PART I done. \n";
}

template <class parameters_type>
void reader<parameters_type>::execute() {
  int nLinesTot = 0;

  std::ifstream input_file(parameters.get_PROCAR_file().c_str());  // input file name

  if (!input_file.is_open())
    std::cout << "\n ERROR: PROCAR_* is not found \n" << std::endl;

  while (std::getline(input_file, oneline))
    ++nLinesTot;

  input_file.clear();

  input_file.seekg(0, input_file.beg);

  int ndp, nstar;

  for (int jj = 0; jj < nLinesTot; jj++) {
    oneline.clear();
    std::getline(input_file, oneline);

    if (!(oneline.empty()) && !(oneline == " ") && !(oneline.find(target1) != std::string::npos) &&
        !(oneline.find(target2) != std::string::npos) &&
        !(oneline.find(target3) != std::string::npos) &&
        !(oneline.find(target4) != std::string::npos) &&
        !(oneline.find(target5) != std::string::npos)) {
      ndp = 0;
      nstar = 0;
      for (int j = 0; j < oneline.length(); j++) {
        if (oneline.substr(j, 1) == ".")
          ndp = ndp + 1;
      }
      for (int j = 0; j < oneline.length(); j++) {
        if (oneline.substr(j, 6) == "******")
          nstar = nstar + 1;
      }

      if (oneline.find(tagk1) != std::string::npos && oneline.find(tagk2) != std::string::npos) {
        oneline.insert(0, "  ");
        oneline.append("  ");
        useful_lines.push_back(oneline);
      }
      else if (oneline.find(tagb1) != std::string::npos && oneline.find(tagb2) != std::string::npos) {
        oneline.insert(0, "  ");
        oneline.append("  ");
        useful_lines.push_back(oneline);
      }
      else if ((ndp + nstar) == 9 || (ndp + nstar) == 10) {
        int ion = 0;

        for (int i = 0; i < oneline.size(); i++)
          if (std::isdigit(oneline[i])) {
            ion = i;
            break;
          }

        if (!(oneline.substr(ion, 1) != "1" || oneline.substr(ion, 1) != "2" ||
              oneline.substr(ion, 1) != "3" || oneline.substr(ion, 1) != "4" ||
              oneline.substr(ion, 1) != "5"))
          std::cout << " Error: atom id problem. " << std::endl;

        if (atof((oneline.substr(ion, 1)).c_str()) == parameters.get_iid()) {
          oneline.insert(0, "  ");
          oneline.append("  ");
          useful_lines.push_back(oneline);
        }
      }
      else {
        std::cout << "\n You missed some possibility. \n" << std::endl;
      }
    }
  }

  std::cout << "\n Total # of useful lines extracted from PROCAR (expected "
            << parameters.get_nKpoints() * (1 + parameters.get_nBands() * 4) << ") is: "
            << " " << useful_lines.size() << std::endl;
  input_file.close();

  int aa;

  // loop k-points
  for (int ik = 0; ik < parameters.get_nKpoints(); ik++) {
    // cut out sub-std::string containing kx,ky,kz
    pos1 = -1;
    pos2 = -1;
    dummy.clear();
    pos1 = useful_lines[ik * (1 + parameters.get_nBands() * 4)].find(":");
    pos2 = useful_lines[ik * (1 + parameters.get_nBands() * 4)].find("weight");
    dummy = useful_lines[ik * (1 + parameters.get_nBands() * 4)].substr(pos1 + 1, pos2 - pos1 - 1);

    // store kx,ky,kz values
    pdp.clear();
    for (int i = 0; i < dummy.length(); i++)
      if (dummy.substr(i, 1) == ".")
        pdp.push_back(i);

    if (pdp.size() != 3)
      std::cout << "\n Error: k-point line does not have 3 d.p. \n" << dummy << std::endl;

    kx(ik) = atof((dummy.substr(pdp[0] - 2, (pdp[1] - 2) - (pdp[0] - 2))).c_str());
    ky(ik) = atof((dummy.substr(pdp[1] - 2, (pdp[2] - 2) - (pdp[1] - 2))).c_str());
    kz(ik) = atof((dummy.substr(pdp[2] - 2)).c_str());

    // std::cout << "\n";
    if (ik < 10)
      std::cout << ik << "\t" << kx(ik) << "\t" << ky(ik) << "\t" << kz(ik) << "\n";
    // std::cout << "\n";

    // loop bands, under that k-point
    for (int ib = 0; ib < parameters.get_nBands(); ib++) {
      int bb = ik * (1 + parameters.get_nBands() * 4) + (1 + ib * 4);

      ndp1 = 0;
      ndp2 = 0;
      ndp3 = 0;
      nstar1 = 0;
      nstar2 = 0;
      nstar3 = 0;
      for (int j = 0; j < useful_lines[bb + 1].length(); j++) {
        if (useful_lines[bb + 1].substr(j, 1) == ".")
          ndp1 = ndp1 + 1;
      }
      for (int j = 0; j < useful_lines[bb + 2].length(); j++) {
        if (useful_lines[bb + 2].substr(j, 1) == ".")
          ndp2 = ndp2 + 1;
      }
      for (int j = 0; j < useful_lines[bb + 3].length(); j++) {
        if (useful_lines[bb + 3].substr(j, 1) == ".")
          ndp3 = ndp3 + 1;
      }
      for (int j = 0; j < useful_lines[bb + 1].length(); j++) {
        if (useful_lines[bb + 1].substr(j, 6) == "******")
          nstar1 = nstar1 + 1;
      }
      for (int j = 0; j < useful_lines[bb + 2].length(); j++) {
        if (useful_lines[bb + 2].substr(j, 6) == "******")
          nstar2 = nstar2 + 1;
      }
      for (int j = 0; j < useful_lines[bb + 3].length(); j++) {
        if (useful_lines[bb + 3].substr(j, 6) == "******")
          nstar3 = nstar3 + 1;
      }

      if (!(useful_lines[bb].find("energy") != std::string::npos) ||
          !(useful_lines[bb].find("occ.") != std::string::npos) || (ndp1 + nstar1) != 10 ||
          (ndp2 + nstar2) != 9 || (ndp3 + nstar3) != 9)
        std::cout << "\n Error: problem in the chunk of (kpoint index, band index)= "
                  << "(" << ik << ", " << ib << ")" << std::endl;

      // deal with band energy line, bb
      pos1 = -1;
      pos2 = -1;
      pos1 = useful_lines[bb].find("energy");
      pos2 = useful_lines[bb].find("occ.");

      band_e(ik, ib) =
          atof((useful_lines[bb].substr(pos1 + 6, (pos2 - 3) - (pos1 + 6) + 1)).c_str()) -
          parameters.get_Fermi_level();  // E_Fermi!
      band_o(ik, ib) = atof((useful_lines[bb].substr(pos2 + 4)).c_str());

      dummy.insert(0, "  ");
      dummy.append("  ");

      {
        // deal with projection magnitude line, bb+1
        aa = bb + 1;
        dummy.clear();
        head.clear();
        tail.clear();
        for (int i = 0; i < (useful_lines[aa].length() - 1); i++)
          if (std::isspace(useful_lines[aa].at(i)) && !std::isspace(useful_lines[aa].at(i + 1)))
            head.push_back(i);
        for (int i = 1; i < (useful_lines[aa].length() - 0); i++)
          if (!std::isspace(useful_lines[aa].at(i - 1)) && std::isspace(useful_lines[aa].at(i)))
            tail.push_back(i);

        if (head.size() != tail.size() || (head.size() != parameters.get_vasp_orbitals() + 2 ||
                                           tail.size() != parameters.get_vasp_orbitals() + 2))  // check
          std::cout << "\n Error: projection magnitude line has problem. " << std::endl;
        if (1.0 * parameters.get_iid() !=
            1.0 * atof((useful_lines[aa].substr(head[0] + 1, tail[0] - head[0] - 1)).c_str()))  // check
          std::cout << "\n Error: you got the wrong impurity atom. " << std::endl;

        for (int j = 0; j < parameters.get_vasp_orbitals(); j++)
          proj_magni(ik, ib, j) =
              1.0 *
              atof((useful_lines[aa].substr(head[j + 1] + 1, tail[j + 1] - head[j + 1] - 1)).c_str());
      }

      {
        // deal with projection phase-Re line, bb+2
        aa = bb + 2;
        dummy.clear();
        head.clear();
        tail.clear();
        for (int i = 0; i < (useful_lines[aa].length() - 1); i++)
          if (std::isspace(useful_lines[aa].at(i)) && !std::isspace(useful_lines[aa].at(i + 1)))
            head.push_back(i);
        for (int i = 1; i < (useful_lines[aa].length() - 0); i++)
          if (!std::isspace(useful_lines[aa].at(i - 1)) && std::isspace(useful_lines[aa].at(i)))
            tail.push_back(i);

        if (head.size() != tail.size() || (head.size() != parameters.get_vasp_orbitals() + 1 ||
                                           tail.size() != parameters.get_vasp_orbitals() + 1))  // check
          std::cout << "\n Error: projection phase-Re line has problem. " << std::endl;
        if (1.0 * parameters.get_iid() !=
            1.0 * atof((useful_lines[aa].substr(head[0] + 1, tail[0] - head[0] - 1)).c_str()))  // check
          std::cout << "\n Error: you got the wrong impurity atom. " << std::endl;

        for (int j = 0; j < parameters.get_vasp_orbitals(); j++)
          proj_phaRe(ik, ib, j) =
              1.0 *
              atof((useful_lines[aa].substr(head[j + 1] + 1, tail[j + 1] - head[j + 1] - 1)).c_str());
      }

      {
        // deal with projection phase-Im line, bb+3
        aa = bb + 3;
        dummy.clear();
        head.clear();
        tail.clear();
        for (int i = 0; i < (useful_lines[aa].length() - 1); i++)
          if (std::isspace(useful_lines[aa].at(i)) && !std::isspace(useful_lines[aa].at(i + 1)))
            head.push_back(i);
        for (int i = 1; i < (useful_lines[aa].length() - 0); i++)
          if (!std::isspace(useful_lines[aa].at(i - 1)) && std::isspace(useful_lines[aa].at(i)))
            tail.push_back(i);

        if (head.size() != tail.size() || (head.size() != parameters.get_vasp_orbitals() + 1 ||
                                           tail.size() != parameters.get_vasp_orbitals() + 1))  // check
          std::cout << "\n Error: projection phase-Im line has problem. " << std::endl;
        if (1.0 * parameters.get_iid() !=
            1.0 * atof((useful_lines[aa].substr(head[0] + 1, tail[0] - head[0] - 1)).c_str()))  // check
          std::cout << "\n Error: you got the wrong impurity atom. " << std::endl;

        for (int j = 0; j < parameters.get_vasp_orbitals(); j++)
          proj_phaIm(ik, ib, j) =
              1.0 *
              atof((useful_lines[aa].substr(head[j + 1] + 1, tail[j + 1] - head[j + 1] - 1)).c_str());
      }
    }  // loop bands

  }  // loop k-points
}

}  // vasp
}  // dft
}  // phys
}  // dca

#endif  // DCA_PHYS_DFT_CONNECTION_VASP_READER_HPP
