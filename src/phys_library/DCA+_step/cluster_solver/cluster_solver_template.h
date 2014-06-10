//-*-C++-*-

#ifndef DCA_CLUSTER_SOLVER_TEMPLATE_H
#define DCA_CLUSTER_SOLVER_TEMPLATE_H

namespace DCA
{
  /*!
   *  \defgroup CLUSTER-SOLVER
   */
  namespace QMCI
  {
    enum QMCI_NAMES {CT_AUX_SOLVER, SS_CT_HYB};

    enum QMCI_SP_MEASUREMENT_NAMES {NFFT};
  }

  /*!
   * \brief   high temperature series expansion solver
   * \author  Peter Staar
   * \version 1.0
   */
  template<CLUSTER_SOLVER_NAMES NAME, LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  class cluster_solver
  {
  public:

    cluster_solver(parameters_type&   parameters_ref,
                   MOMS_type&         MOMS_ref);

    ~cluster_solver();

    void initialize(int dca_iteration);

    void execute();

    template<typename dca_info_struct_t>
    void finalize(dca_info_struct_t& dca_info_struct);


    void read(std::string filename);

    void write(std::string filename);

    template<IO::FORMAT DATA_FORMAT>
    void read(IO::reader<DATA_FORMAT>& reader);

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);
  };

}

#endif
