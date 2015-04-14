//-*-C++-*-

#ifndef MATERIAL_INTERACTION_TEMPLATE_H
#define MATERIAL_INTERACTION_TEMPLATE_H

/*!
 *  \author Peter Staar
 */
template<material_name_type name, typename point_group_type>
class material_interaction
{
public:

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
                                       parameters_type&            parameters);
};

#endif
