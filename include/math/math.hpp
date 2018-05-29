//#pragma once

#ifndef __MATHJOY_MATH_HPP__
#define __MATHJOY_MATH_HPP__

template <typename N>
inline bool odd(N n) {
  return (n & 0x1);
}

/*
 * f could be the operation representing addition or multiplication
 */
template <typename R, typename N, typename BinaryOperation>
R power_recursive1(R acc, R Base, N n, BinaryOperation f) {
  // Haven't thought of a good way to provide an identity element
  if ( n < 1 ) throw;

  if ( n == 1 ) return acc;
  if ( odd(n) ) {
    return power_recursive(f(acc, Base), f(Base, Base), n >> 1, f);
  }
  else {
    return power_recursive(acc, f(Base, Base), n >> 1, f);
  }
}

/*
 * Strict tail recursion
 */
template <typename R, typename N, typename BinaryOperation>
R power_recursive2(R acc, R Base, N n, BinaryOperation f) {
  // Haven't thought of a good way to provide an identity element
  if ( n < 1 ) throw;

  if ( odd(n) ) {
    acc = f(acc, Base);
    if ( n == 1 ) return acc;
  }
  Base = f(Base, Base);
  n >>= 1;
  return power_recursive(acc, Base, n, f);
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
