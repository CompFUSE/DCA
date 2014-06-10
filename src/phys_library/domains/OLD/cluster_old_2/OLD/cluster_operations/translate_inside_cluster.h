//-*-C++-*-

#ifndef TRANSLATE_INSIDE_CLUSTER_H
#define TRANSLATE_INSIDE_CLUSTER_H

namespace CLUSTER_OPERATIONS
{
  /*!
   *  \author Peter Staar
   */
  class translate_inside_cluster
  {
  public:
    
    template<typename scalar_type>
    static std::vector<scalar_type> execute(std::vector<scalar_type>& r,
					    std::vector<std::vector<scalar_type> >& B);
  };
  
  template<typename scalar_type>
  std::vector<scalar_type> translate_inside_cluster::execute(std::vector<scalar_type>& r,
							     std::vector<std::vector<scalar_type> >& B)
  {
    int DIMENSION = r.size();

    std::vector<double> r_affine = VECTOR_OPERATIONS::COORDINATES(r, B);

    for(size_t d=0; d<r.size(); d++)
      {
	while(r_affine[d]<-1.e-6)
	  r_affine[d] += 1.;
    
	while(r_affine[d]>1-1.e-6)
	  r_affine[d] -= 1.;
      }

    std::vector<double> r_vec(r.size(), 0.);
    
    for(size_t d1=0; d1<DIMENSION; ++d1)
      for(size_t d0=0; d0<DIMENSION; ++d0)
	r_vec[d0] += r_basis[d1][d0]*r_affine[d1];
    
    return r_vec;
  }

}

#endif
