//-*-C++-*-

/*
 *      Author: Peter Staar
 */


#ifndef D4_3D_H
#define D4_3D_H

class D4_3D : public point_group<3>
{
public:

  const static int size = 9;

  typedef P_3D inversion_type;

  typedef Cn_3D<0,0,1,1,4> Cn_3D_Z_1_4_type;
  typedef Cn_3D<0,0,1,2,4> Cn_3D_Z_2_4_type;
  typedef Cn_3D<0,0,1,3,4> Cn_3D_Z_3_4_type;

  typedef Sn_3D<2,0,8> Sn_3D_0_8_type;
  typedef Sn_3D<2,1,8> Sn_3D_1_8_type;
  typedef Sn_3D<2,2,8> Sn_3D_2_8_type;
  typedef Sn_3D<2,3,8> Sn_3D_3_8_type;
  typedef Sn_3D<2,4,8> Sn_3D_4_8_type;

  typedef Typelist<inversion_type,
		      
		      Cn_3D_Z_1_4_type,
		      Cn_3D_Z_2_4_type,
		      Cn_3D_Z_3_4_type,
		      
		      Sn_3D_0_8_type,
		      Sn_3D_1_8_type,
		      Sn_3D_2_8_type,
		      Sn_3D_3_8_type,
		      Sn_3D_4_8_type> point_group_type_list;
};


#endif
