//-*-C++-*-

#ifndef GENERIC_ASSERT_H
#define GENERIC_ASSERT_H

using namespace TL;  

/*!
 *  \author: Peter Staar
 */
template<bool>
struct GENERIC_ASSERT
{
  static void execute(){}
};

template<>
struct GENERIC_ASSERT<false>
{};

#endif
