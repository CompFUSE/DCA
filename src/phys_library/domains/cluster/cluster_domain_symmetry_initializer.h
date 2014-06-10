//-*-C++-*-

#ifndef INITIALIZE_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_H
#define INITIALIZE_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_H

template<typename cluster_type, typename point_group_type>
class cluster_domain_symmetry_initializer
{
  typedef typename cluster_symmetry<cluster_type>::cluster_family_type cluster_family_type;

public:

  static void execute()
  {
    cluster_reduction<cluster_family_type, point_group_type> cluster_reduction_obj;
    cluster_reduction_obj.execute();
  }
};

template<typename cluster_type, typename point_group_type>
class cluster_domain_symmetry_initializer<dmn_0<cluster_type>, point_group_type>
{
  typedef typename cluster_symmetry<cluster_type>::cluster_family_type cluster_family_type;

public:

  static void execute()
  {
    cluster_reduction<cluster_family_type, point_group_type> cluster_reduction_obj;
    cluster_reduction_obj.execute();
  }
};


#endif 
