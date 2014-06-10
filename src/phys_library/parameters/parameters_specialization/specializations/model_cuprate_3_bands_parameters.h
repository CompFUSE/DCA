// //-*-C++-*-

// /*
//  *      Author: peter staar
//  */


// #ifndef MODEL_CUPRATE_3_BANDS_PARAMETERS_H
// #define MODEL_CUPRATE_3_BANDS_PARAMETERS_H

// template<typename DCA_point_group_type, typename interaction_t>
// class model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >
// {
// public:

//   model_parameters();
//   ~model_parameters();

// /******************************************
//  ***        CONCURRENCY                 ***
//  ******************************************/

//   template<class concurrency_type>
//   int  get_buffer_size(const concurrency_type& concurrency) const;

//   template<class concurrency_type>
//   void pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

//   template<class concurrency_type>
//   void unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

// /******************************************
//  ***        READ/WRITE                  ***
//  ******************************************/

//   template<class stream_type>
//   void to_JSON(stream_type& ss, bool is_end=false);
  
//   template<class JSON_reader_type>
//   void from_JSON(JSON_reader_type& reader);

// /******************************************
//  ***        DATA                        ***
//  ******************************************/

//   double  get_U_hubbard();
//   double  get_J_exchange();

//   std::vector<int>& get_H_k_grid_size();
  
//   std::vector<complex<double> > get_dd();
//   std::vector<complex<double> > get_pd();
//   std::vector<complex<double> > get_pp();

//  private:

//   double U_hubbard;
//   double J_exchange;

//   std::vector<int> H_k_grid_size;

//   std::vector<complex<double> > dd;
//   std::vector<complex<double> > pd;
//   std::vector<complex<double> > pp;
// };

// template<typename DCA_point_group_type, typename interaction_t>
// model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::model_parameters()
// {}

// template<typename DCA_point_group_type, typename interaction_t>
// model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::~model_parameters()
// {}

// /******************************************
//  ***        CONCURRENCY                 ***
//  ******************************************/

// template<typename DCA_point_group_type, typename interaction_t>
// template<class concurrency_type>
// int model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_buffer_size(const concurrency_type& concurrency) const
// {
//   int buffer_size = 0;

//   buffer_size += concurrency.getBufferSize(U_hubbard);
//   buffer_size += concurrency.getBufferSize(J_exchange);
//   buffer_size += concurrency.getBufferSize(H_k_grid_size);

//   buffer_size += concurrency.getBufferSize(dd);
//   buffer_size += concurrency.getBufferSize(pd);
//   buffer_size += concurrency.getBufferSize(pp);

//   return buffer_size;
// }

// template<typename DCA_point_group_type, typename interaction_t>
// template<class concurrency_type>
// void model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
// {
//   concurrency.pack(buffer, buffer_size, position, U_hubbard);
//   concurrency.pack(buffer, buffer_size, position, J_exchange);
//   concurrency.pack(buffer, buffer_size, position, H_k_grid_size);

//   concurrency.pack(buffer, buffer_size, position, dd);
//   concurrency.pack(buffer, buffer_size, position, pd);
//   concurrency.pack(buffer, buffer_size, position, pp);
// }

// template<typename DCA_point_group_type, typename interaction_t>
// template<class concurrency_type>
// void model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
// {
//   concurrency.unpack(buffer, buffer_size, position, U_hubbard);
//   concurrency.unpack(buffer, buffer_size, position, J_exchange);
//   concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);

//   concurrency.unpack(buffer, buffer_size, position, dd);
//   concurrency.unpack(buffer, buffer_size, position, pd);
//   concurrency.unpack(buffer, buffer_size, position, pp);
// }

// /******************************************
//  ***        READ/WRITE                  ***
//  ******************************************/

// template<typename DCA_point_group_type, typename interaction_t>
// template<class stream_type>
// void model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::to_JSON(stream_type& ss, bool is_end)
// {
//   ss << "\"model-parameters\" :";
//   ss << "\n{ \n";

//   JSON_writer::write(ss, "U_hubbard"                , U_hubbard);
//   JSON_writer::write(ss, "J_exchange"               , J_exchange);
  
//   JSON_writer::write(ss, "dd", dd);
//   JSON_writer::write(ss, "pd", pd);
//   JSON_writer::write(ss, "pp", pp);

//   JSON_writer::write(ss, "H_k_grid_size"            , H_k_grid_size, true);

//   if(is_end)
//     ss << "}\n";
//   else
//     ss << "},\n";
// }
 
// template<typename DCA_point_group_type, typename interaction_t> 
// template<class JSON_reader_type>
// void model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::from_JSON(JSON_reader_type& reader)
// {
//   typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
//   const JsonAccessor control(reader["model-parameters"]);

//   U_hubbard          <= control["U_hubbard"];
//   J_exchange         <= control["J_exchange"];
//   H_k_grid_size      <= control["H_k_grid_size"]; 

//   dd <= control["dd"];
//   pd <= control["pd"];
//   pp <= control["pp"];

// //   for(size_t l=0; l<dd.size(); l++)
// //     cout << dd[l];
// //   cout << endl;
// //   for(size_t l=0; l<pd.size(); l++)
// //     cout << pd[l];
// //   cout << endl;
// //   for(size_t l=0; l<pp.size(); l++)
// //     cout << pp[l];
// //   cout << endl;
// //   throw std::logic_error("STOP");

//   if(int(H_k_grid_size.size()) != tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t>::DIMENSION)  
//     throw std::logic_error("int(H_k_grid_size.size()) != model::DIMENSION");
// }

// /******************************************
//  ***        DATA                        ***
//  ******************************************/

// template<typename DCA_point_group_type, typename interaction_t> 
// double model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_U_hubbard()
// {
//   return U_hubbard;
// }
 
// template<typename DCA_point_group_type, typename interaction_t> 
// double model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_J_exchange()
// {
//   return J_exchange;
// }

// template<typename DCA_point_group_type, typename interaction_t> 
// std::vector<int>& model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_H_k_grid_size()
// {
//   return H_k_grid_size;
// }

// template<typename DCA_point_group_type, typename interaction_t> 
// std::vector<complex<double> > model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_dd()
// {
//   return dd;
// }
 
// template<typename DCA_point_group_type, typename interaction_t> 
// std::vector<complex<double> > model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_pd()
// {
//   return pd;
// }

// template<typename DCA_point_group_type, typename interaction_t> 
// std::vector<complex<double> > model_parameters<tight_binding_model<cuprate_model_3_bands<DCA_point_group_type>, interaction_t> >::get_pp()
// {
//   return pp;
// }
// #endif
