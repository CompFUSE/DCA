//-*-C++-*-

/*
 *      Author: peterstaar
 */


#ifndef ROT_4_H_
#define ROT_4_H_

class Rot_4 : public point_group<2>
{
public:

  const static int size = 3;

  typedef Cn_2D<1,4> Cn_2D_1_4_type;
  typedef Cn_2D<2,4> Cn_2D_2_4_type;
  typedef Cn_2D<3,4> Cn_2D_3_4_type;

  typedef TYPELIST_3(Cn_2D_1_4_type, 
		     Cn_2D_2_4_type, 
		     Cn_2D_3_4_type) symmetry_list;
};


#endif
