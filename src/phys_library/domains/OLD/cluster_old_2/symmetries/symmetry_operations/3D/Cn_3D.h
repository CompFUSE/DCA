//-*-C++-*-

#ifndef Cn_3D_H_
#define Cn_3D_H_

/*!
 *  \Author Peter Staar
 *  \brief  Elements of the cyclic group C_n(m)
 */
template<int ux, int uy, int uz, int n, int m>
class Cn_3D : public group_action<3>
{
public:

  typedef group_action<3>  base_type;
  typedef Cn_3D< ux, uy, uz,n,m>  this_type;

  Cn_3D(){};
  ~Cn_3D(){};

  static double* matrix()
  {
    // rotation around an axis {ux,uy,uz} with angle th = 2*pi*n/m;

    double c = COSINE_EVAL<n,m>::value();//cos(2*M_PI*double(n)/double(m));
    double s = SINE_EVAL<n,m>::value();//sin(2*M_PI*double(n)/double(m));

    double Ux = double(ux)/sqrt(double(ux*ux+uy*uy+uz*uz));
    double Uy = double(uy)/sqrt(double(ux*ux+uy*uy+uz*uz));
    double Uz = double(uz)/sqrt(double(ux*ux+uy*uy+uz*uz));
    
    static double matrix[3*3] = { Ux*Ux + (1-Ux*Ux)*c, Ux*Uy*(1-c)+Uz*s    , Ux*Uz*(1-c)-Uy*s,
				  Ux*Uy*(1-c)-Uz*s   , Uy*Uy + (1-Uy*Uy)*c , Uy*Uz*(1-c)+Ux*s,
				  Ux*Uz*(1-c)+Uy*s   , Uy*Uz*(1-c)-Ux*s    , Uz*Uz+(1-Uz*Uz)*c};
    return matrix;
  }

private:
  
  
};

#endif
