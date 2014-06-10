//-*-C++-*-


// elements of the cyclic group C_n(m)

#ifndef MIRR_ALPHA_BETA_3D_H_
#define MIRR_ALPHA_BETA_3D_H_

/*!
 *   \author Peter Staar
 */
template<int alpha_num, int alpha_den,
	 int beta_num , int beta_den>
class mirr_a_b : public group_action<3>
{
public:

  typedef group_action<3>  base_type;
  typedef mirr_a_b< alpha_num, alpha_den,
		    beta_num , beta_den> this_type;

  mirr_a_b(){};
  ~mirr_a_b(){};

  const static double* matrix()
  {
     double a = 2*M_PI*double(alpha_num)/double(alpha_den); // polar coordinate \theta
     double b = 2*M_PI*double(beta_num) /double(beta_den);  // polar coordinate \phi 
     
    // matrix in column-major format...
     static double matrix[3*3] = {pow(cos(a),2) - sin(a)*(-(pow(cos(b),2)*sin(a)) + sin(a)*pow(sin(b),2)),
				       cos(a)*sin(a) - sin(a)*(cos(a)*pow(cos(b),2) - cos(a)*pow(sin(b),2)),
				       -2*cos(b)*sin(a)*sin(b),

				       cos(a)*sin(a) + cos(a)*(-(pow(cos(b),2)*sin(a)) + sin(a)*pow(sin(b),2)),
				       pow(sin(a),2) + cos(a)*(cos(a)*pow(cos(b),2) - cos(a)*pow(sin(b),2)),
				       2*cos(a)*cos(b)*sin(b),

				       -2*cos(b)*sin(a)*sin(b),
				       2*cos(a)*cos(b)*sin(b),
				       -pow(cos(b),2) + pow(sin(b),2)};
    return matrix;
  }

  /*
  // MATHEMATICA CODE

  // Z[\[Alpha]_] := ({
  //    {Cos[\[Alpha]], -Sin[\[Alpha]], 0},
  //    {Sin[\[Alpha]], Cos[\[Alpha]], 0},
  //    {0, 0, 1}
  //   })
  // X[\[Alpha]_] := ({
  //    {1, 0, 0},
  //    {0, Cos[\[Alpha]], -Sin[\[Alpha]]},
  //    {0, Sin[\[Alpha]], Cos[\[Alpha]]}
  //   })
  // P := ({
  //    {1, 0, 0},
  //    {0, 1, 0},
  //    {0, 0, -1}
  //   })
  //
  // Transpose[
  //   Z[\[Alpha]].X[\[Beta]].P.X[-\[Beta]].Z[-\[Alpha]]] 
  // matrix = CForm[Transpose[
  //    X[\[Alpha]].Z[\[Beta]].P.Z[-\[Beta]].X[-\[Alpha]]] /. {\
  // \[Alpha] -> a, \[Beta] -> b, }]
  */

};

#endif
