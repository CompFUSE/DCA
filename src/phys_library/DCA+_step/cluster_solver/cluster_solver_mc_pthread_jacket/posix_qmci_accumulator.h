//-*-C++-*-

#ifndef DCA_QMCI_POSIX_JACKET_FOR_MC_ACCUMULATION_H
#define DCA_QMCI_POSIX_JACKET_FOR_MC_ACCUMULATION_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     * \ingroup POSIX-TEMPLATES
     * \brief   A posix-jacket that implements a MC-accumulator, independent of the MC-TYPE.
     * \author  Peter Staar, Raffaele Solca
     * \version 1.0
     */
    template<class qmci_accumulator_type>
    class posix_qmci_accumulator : protected qmci_accumulator_type
    {
      typedef typename qmci_accumulator_type::my_parameters_type   parameters_type;
      typedef typename qmci_accumulator_type::my_MOMS_type         MOMS_type;

      typedef posix_qmci_accumulator<qmci_accumulator_type>          this_type;

    public:

      posix_qmci_accumulator(parameters_type& parameters_ref,
                             MOMS_type&       MOMS_ref,
                             int              id);

      ~posix_qmci_accumulator();

      using qmci_accumulator_type::initialize;
      using qmci_accumulator_type::finalize;
      //using qmci_accumulator_type::to_JSON;
      using qmci_accumulator_type::get_configuration;

      template<typename walker_type>
      void update_from(walker_type& walker);

      void wait_for_qmci_walker();

      void measure(pthread_mutex_t&        mutex_queue,
                   std::queue<this_type*>& accumulators_queue);

      //void sum_to(qmci_accumulator_type& accumulator_obj);
      //int get_expansion_order();

      template<LIN_ALG::device_type device_t, class parameters_t, class MOMS_t>
      void sum_to(MC_accumulator<CT_AUX_SOLVER, device_t, parameters_t, MOMS_t>& accumulator_obj);

      template<LIN_ALG::device_type device_t, class parameters_t, class MOMS_t>
      void sum_to(MC_accumulator<SS_CT_HYB, device_t, parameters_t, MOMS_t>& accumulator_obj);

    protected:

      using qmci_accumulator_type::get_Gflop;
      using qmci_accumulator_type::get_sign;
      using qmci_accumulator_type::get_number_of_measurements;

    private:

      using qmci_accumulator_type::parameters;
      using qmci_accumulator_type::MOMS;


      int             thread_id;

      bool            measuring;
      pthread_cond_t  start_measuring;

      pthread_mutex_t mutex_accumulator;
    };

    template<class qmci_accumulator_type>
    posix_qmci_accumulator<qmci_accumulator_type>::posix_qmci_accumulator(parameters_type& parameters_ref,
                                                                          MOMS_type&       MOMS_ref,
                                                                          int              id):
      qmci_accumulator_type(parameters_ref, MOMS_ref, id),
      thread_id(id),
      measuring(false)
    {
      pthread_cond_init (&start_measuring  , NULL);

      pthread_mutex_init(&mutex_accumulator      , NULL);
    }

    template<class qmci_accumulator_type>
    posix_qmci_accumulator<qmci_accumulator_type>::~posix_qmci_accumulator()
    {
      pthread_cond_destroy (&start_measuring);

      pthread_mutex_destroy(&mutex_accumulator);
    }


    template<class qmci_accumulator_type>
    template<typename walker_type>
    void posix_qmci_accumulator<qmci_accumulator_type>::update_from(walker_type& walker)
    {
      {
        pthread_mutex_lock(&mutex_accumulator);

        if(measuring) throw std::logic_error(__FUNCTION__);

        qmci_accumulator_type::update_from(walker);

        measuring = true;

        pthread_cond_signal(&start_measuring);

        pthread_mutex_unlock(&mutex_accumulator);
      }
    }

    template<class qmci_accumulator_type>
    void posix_qmci_accumulator<qmci_accumulator_type>::wait_for_qmci_walker()
    {
      pthread_mutex_lock(&mutex_accumulator);

      while(!measuring)
        pthread_cond_wait(&start_measuring, &mutex_accumulator);

      pthread_mutex_unlock(&mutex_accumulator);
    }

    template<class qmci_accumulator_type>
    void posix_qmci_accumulator<qmci_accumulator_type>::measure(pthread_mutex_t&        /*mutex_queue*/,
                                                                std::queue<this_type*>& /*accumulators_queue*/)
    {
      pthread_mutex_lock(&mutex_accumulator);

      qmci_accumulator_type::measure();

      measuring = false;

      pthread_mutex_unlock(&mutex_accumulator);
    }

    template<class qmci_accumulator_type>
    template<LIN_ALG::device_type device_t, class parameters_t, class MOMS_t>
    void posix_qmci_accumulator<qmci_accumulator_type>::sum_to(MC_accumulator<CT_AUX_SOLVER, device_t, parameters_t, MOMS_t>& accumulator_obj)
    {
      typedef MC_accumulator<CT_AUX_SOLVER, device_t, parameters_t, MOMS_t> current_qmci_accumulator_type;

      pthread_mutex_lock(&mutex_accumulator);

      finalize();

      accumulator_obj.get_Gflop()                  += current_qmci_accumulator_type::get_Gflop();

      accumulator_obj.get_sign()                   += current_qmci_accumulator_type::get_sign();
      accumulator_obj.get_number_of_measurements() += current_qmci_accumulator_type::get_number_of_measurements();

      {
	for(int i=0; i<current_qmci_accumulator_type::visited_expansion_order_k.size(); i++)
	  accumulator_obj.get_visited_expansion_order_k()(i) += current_qmci_accumulator_type::visited_expansion_order_k(i);

	for(int i=0; i<current_qmci_accumulator_type::error.size(); i++)
	  accumulator_obj.get_error_distribution()(i) += current_qmci_accumulator_type::error(i);
      }

      {// equal time measurements
	for(int i=0; i<current_qmci_accumulator_type::G_r_t.size(); i++)
	  accumulator_obj.get_G_r_t()(i) += current_qmci_accumulator_type::G_r_t(i);
	
	for(int i=0; i<current_qmci_accumulator_type::G_r_t_stddev.size(); i++)
	  accumulator_obj.get_G_r_t_stddev()(i) += current_qmci_accumulator_type::G_r_t_stddev(i);
	
	for(int i=0; i<current_qmci_accumulator_type::charge_cluster_moment.size(); i++)
	  accumulator_obj.get_charge_cluster_moment()(i) += current_qmci_accumulator_type::charge_cluster_moment(i);
	
	for(int i=0; i<current_qmci_accumulator_type::magnetic_cluster_moment.size(); i++)
	  accumulator_obj.get_magnetic_cluster_moment()(i) += current_qmci_accumulator_type::magnetic_cluster_moment(i);
	
	for(int i=0; i<current_qmci_accumulator_type::dwave_pp_correlator.size(); i++)
	  accumulator_obj.get_dwave_pp_correlator()(i) += current_qmci_accumulator_type::dwave_pp_correlator(i);
      }

      {// sp-measurements
	for(int i=0; i<current_qmci_accumulator_type::K_r_t.size(); i++)
	  accumulator_obj.get_K_r_t()(i) += current_qmci_accumulator_type::K_r_t(i);

	for(int i=0; i<current_qmci_accumulator_type::M_r_w.size(); i++)
	  accumulator_obj.get_M_r_w()(i) += current_qmci_accumulator_type::M_r_w(i);
	
	for(int i=0; i<current_qmci_accumulator_type::M_r_w_squared.size(); i++)
	  accumulator_obj.get_M_r_w_squared()(i) += current_qmci_accumulator_type::M_r_w_squared(i);
      }

      {// tp-measurements
	for(int i=0; i<current_qmci_accumulator_type::G4.size(); i++)
	  accumulator_obj.get_G4()(i) += current_qmci_accumulator_type::G4(i);
      }

      pthread_mutex_unlock(&mutex_accumulator);
    }

    template<class qmci_accumulator_type>
    template<LIN_ALG::device_type device_t, class parameters_t, class MOMS_t>
    void posix_qmci_accumulator<qmci_accumulator_type>::sum_to(MC_accumulator<SS_CT_HYB, device_t, parameters_t, MOMS_t>& accumulator_obj)
    {
      typedef MC_accumulator<SS_CT_HYB, device_t, parameters_t, MOMS_t> current_qmci_accumulator_type;

      pthread_mutex_lock(&mutex_accumulator);

      finalize();

      accumulator_obj.get_sign()                   += current_qmci_accumulator_type::get_sign();
      accumulator_obj.get_number_of_measurements() += current_qmci_accumulator_type::get_number_of_measurements();

      for(int i=0; i<current_qmci_accumulator_type::visited_expansion_order_k.size(); i++)
        accumulator_obj.get_visited_expansion_order_k()(i) += current_qmci_accumulator_type::visited_expansion_order_k(i);

      {// sp-measurements
	for(int i=0; i<current_qmci_accumulator_type::G_r_w.size(); i++)
	  accumulator_obj.get_G_r_w()(i) += current_qmci_accumulator_type::G_r_w(i);

	for(int i=0; i<current_qmci_accumulator_type::GS_r_w.size(); i++)
	  accumulator_obj.get_GS_r_w()(i) += current_qmci_accumulator_type::GS_r_w(i);
      }

      pthread_mutex_unlock(&mutex_accumulator);
    }


  }//namespace QMCI

}//namespace DCA

#endif
