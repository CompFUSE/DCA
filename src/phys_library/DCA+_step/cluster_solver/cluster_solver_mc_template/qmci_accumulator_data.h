//-*-C++-*-

#ifndef DCA_QMCI_ACCUMULATOR_DATA_H
#define DCA_QMCI_ACCUMULATOR_DATA_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     * \class   MC_bare_data
     * \brief   Base class containing bare data for a MC accumulator;
     * \author  Peter Staar
     * \version 1.0
     */
    class MC_accumulator_data
    {
    public:

      double& get_Gflop();

      double& get_sign();
      double& get_number_of_measurements();

      void initialize(int dca_iteration);

    protected:

      int DCA_iteration;

      double GFLOP;

      double current_sign;
      double accumulated_sign;

      double number_of_measurements;
    };

    double& MC_accumulator_data::get_Gflop()
    {
      return GFLOP;
    }

    double& MC_accumulator_data::get_sign()
    {
      return accumulated_sign;
    }

    double& MC_accumulator_data::get_number_of_measurements()
    {
      return number_of_measurements;
    }

    void MC_accumulator_data::initialize(int dca_iteration)
    {
      DCA_iteration = dca_iteration;

      GFLOP = 0.;

      current_sign     = 1;
      accumulated_sign = 0;

      number_of_measurements = 0;
    }

  }

}

#endif
