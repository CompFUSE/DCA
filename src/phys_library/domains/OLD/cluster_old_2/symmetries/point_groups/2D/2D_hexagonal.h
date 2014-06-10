//-*-C++-*-

/*
 *      Author: Peter Staar
 */


#ifndef HEXAGONAL_2D_H_
#define HEXAGONAL_2D_H_

/*
 *  group-actions
 */

typedef Cn_2D<1,3> Cn_2D_1_3_type;
typedef Cn_2D<2,3> Cn_2D_2_3_type;

typedef Sn_2D<0,3> Sn_2D_0_3_type;
typedef Sn_2D<1,3> Sn_2D_1_3_type;
typedef Sn_2D<2,3> Sn_2D_2_3_type;

typedef Cn_2D<1,6> Cn_2D_1_6_type;
typedef Cn_2D<2,6> Cn_2D_2_6_type;
typedef Cn_2D<3,6> Cn_2D_3_6_type;
typedef Cn_2D<4,6> Cn_2D_4_6_type;
typedef Cn_2D<5,6> Cn_2D_5_6_type;

typedef Sn_2D<0,6> Sn_2D_0_6_type;
typedef Sn_2D<1,6> Sn_2D_1_6_type;
typedef Sn_2D<2,6> Sn_2D_2_6_type;
typedef Sn_2D<3,6> Sn_2D_3_6_type;
typedef Sn_2D<4,6> Sn_2D_4_6_type;
typedef Sn_2D<5,6> Sn_2D_5_6_type;

/*
 *  pointgroup :: set of group-actions
 */

typedef TYPELIST_2(Cn_2D_1_3_type, 
		   Cn_2D_2_3_type) C3;

typedef TYPELIST_3(Sn_2D_0_3_type,
		   Sn_2D_1_3_type, 
		   Sn_2D_2_3_type) S3;

typedef TYPELIST_5(Cn_2D_1_3_type, 
		   Cn_2D_2_3_type, 
		   Sn_2D_0_3_type,
		   Sn_2D_1_3_type, 
		   Sn_2D_2_3_type) D3;

typedef TYPELIST_5(Cn_2D_1_6_type, 
		   Cn_2D_2_6_type, 
		   Cn_2D_3_6_type,
		   Cn_2D_4_6_type,
		   Cn_2D_5_6_type) C6;

typedef TYPELIST_6(Sn_2D_0_6_type,
		   Sn_2D_1_6_type, 
		   Sn_2D_2_6_type, 
		   Sn_2D_3_6_type, 
		   Sn_2D_4_6_type,
		   Sn_2D_5_6_type) S6;

typedef TYPELIST_11(Cn_2D_1_6_type, 
		    Cn_2D_2_6_type, 
		    Cn_2D_3_6_type,
		    Cn_2D_4_6_type,
		    Cn_2D_5_6_type,
		    Sn_2D_0_6_type,
		    Sn_2D_1_6_type, 
		    Sn_2D_2_6_type, 
		    Sn_2D_3_6_type, 
		    Sn_2D_4_6_type,
		    Sn_2D_5_6_type) D6;

#endif
