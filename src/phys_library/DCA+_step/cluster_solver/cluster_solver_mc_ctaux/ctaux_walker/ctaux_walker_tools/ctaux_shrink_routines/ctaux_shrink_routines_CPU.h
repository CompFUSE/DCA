//-*-C++-*-

#ifndef DCA_QMCI_SHRINK_ROUTINES_CPU_H
#define DCA_QMCI_SHRINK_ROUTINES_CPU_H

namespace DCA
{
  namespace QMCI
  {
    template<>
    class SHRINK_TOOLS_ALGORITHMS<LIN_ALG::CPU>
    {
    public:

      SHRINK_TOOLS_ALGORITHMS(int id):
        thread_id(id),
        stream_id(0)
      {}

      ~SHRINK_TOOLS_ALGORITHMS()
      {}

      void execute(std::vector<int>& source_index,
                   std::vector<int>& target_index,
                   LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
                   LIN_ALG::matrix<double, LIN_ALG::CPU>& G0)
      {
        assert(source_index.size()==target_index.size());

        assert(N.get_number_of_rows() == G0.get_number_of_rows());
        assert(N.get_number_of_cols() == G0.get_number_of_cols());

        for(size_t l=0; l<source_index.size(); ++l){
          LIN_ALG::SWAP<LIN_ALG::CPU>::row_and_column(N , source_index[l], target_index[l], thread_id, stream_id);
          LIN_ALG::SWAP<LIN_ALG::CPU>::row_and_column(G0, source_index[l], target_index[l], thread_id, stream_id);
        }
      }

      void execute(std::vector<int>& source_index_up,
                   std::vector<int>& target_index_up,
                   LIN_ALG::matrix<double, LIN_ALG::CPU>& N_up,
                   LIN_ALG::matrix<double, LIN_ALG::CPU>& G0_up,
                   std::vector<int>& source_index_dn,
                   std::vector<int>& target_index_dn,
                   LIN_ALG::matrix<double, LIN_ALG::CPU>& N_dn,
                   LIN_ALG::matrix<double, LIN_ALG::CPU>& G0_dn)
      {
        assert(source_index_up.size()==target_index_up.size());
        assert(source_index_dn.size()==target_index_dn.size());

        assert(N_up.get_current_size() == G0_up.get_current_size());
        assert(N_dn.get_current_size() == G0_dn.get_current_size());

        for(size_t l=0; l<source_index_up.size(); ++l)
          LIN_ALG::SWAP<LIN_ALG::CPU>::row_and_column(N_up, source_index_up[l], target_index_up[l], thread_id, stream_id);

        for(size_t l=0; l<source_index_up.size(); ++l)
          LIN_ALG::SWAP<LIN_ALG::CPU>::row_and_column(G0_up, source_index_up[l], target_index_up[l], thread_id, stream_id);

        for(size_t l=0; l<source_index_dn.size(); ++l)
          LIN_ALG::SWAP<LIN_ALG::CPU>::row_and_column(N_dn, source_index_dn[l], target_index_dn[l], thread_id, stream_id);

        for(size_t l=0; l<source_index_dn.size(); ++l)
          LIN_ALG::SWAP<LIN_ALG::CPU>::row_and_column(G0_dn, source_index_dn[l], target_index_dn[l], thread_id, stream_id);
      }

    private:

      int thread_id;
      int stream_id;
    };

  }

}

#endif
