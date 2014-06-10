//-*-C++-*-

#ifndef IS_EQUAL_TYPE_H
#define IS_EQUAL_TYPE_H

using namespace std; 
using namespace TL;  

/*!
 *  \author Peter Staar
 */
template<typename type_1, typename type_2>
struct IS_EQUAL_TYPE
{
  const static bool check = false;
};

template<typename type_1>
struct IS_EQUAL_TYPE<type_1, type_1>
{
  const static bool check = true;
};

#endif
