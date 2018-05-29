//#pragma once

#ifndef __MATHJOY_MATH_HPP__
#define __MATHJOY_MATH_HPP__

template <typename N>
inline bool odd(N n) {
  return (n & 0x1);
}
/*
 * A generic power function
 * Only accept n >= 1
 */
template <typename R, typename N, typename BinaryOperation>
R power(R Base, N n, BinaryOperation f) {
  if ( n < 1 ) throw;

  R acc = Base;

  while ( true ) {
    if ( odd(n) ) {
      acc = f(acc, Base);
      if ( n == 1 )
        return acc;
    }

    // half n
    n >>= 1;
    Base = fcc(Base, Base);
  }
}

#endif
