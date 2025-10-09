// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// ADIOS2 global object


/** For testing only, in main_dca adios2::ADIOS is a member of concurrency the concurrency context owned by main.
 */
class GlobalAdios
{
private:
  GlobalAdios();
public:
  static adios2::ADIOS& getAdios()
  {
    static adios2::ADIOS adios_("", MPI_COMM_SELF);
    return adios_;
  }
};
