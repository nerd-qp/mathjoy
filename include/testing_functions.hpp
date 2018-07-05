
#ifndef _TESTING_FUNCTIONS_
#define _TESTING_FUNCTIONS_


#include <cmath>
#define WHICH_FUNC 2

double test_p(double x) {
#if WHICH_FUNC == 1
  return 1;
#elif WHICH_FUNC == 2
  return 1;
#elif WHICH_FUNC == 3
  return 1;
#elif WHICH_FUNC == 4
  return 0;
#elif WHICH_FUNC == 5
  return x*x/11;
#else
  return -1;
#endif
}

double test_f(double x) {
#if WHICH_FUNC == 1
  return -1;
#elif WHICH_FUNC == 2
  return M_PI * M_PI / 2 * sin(M_PI * x);
#elif WHICH_FUNC == 3
  return x*x - 2*x - 2;
#elif WHICH_FUNC == 4
  return x;
#elif WHICH_FUNC == 5
  return - 2 * 1e3 / 11 * pow(x, 9);
#else
  return -x;
#endif
}

double test_q(double x) {
#if WHICH_FUNC == 1
  return 0;
#elif WHICH_FUNC == 2
  return M_PI * M_PI / 4;
#elif WHICH_FUNC == 3
  return 1;
#elif WHICH_FUNC == 4
  return 1;
#elif WHICH_FUNC == 5
  return 10;
#else
  return 1;
#endif
}

double result_f(double x) {
#if WHICH_FUNC == 1
  return x*x/2 - x;
#elif WHICH_FUNC == 2
  return sin(M_PI/2 * x);
#elif WHICH_FUNC == 3
  return x*x-2*x;
#elif WHICH_FUNC == 4
  return x;
#elif WHICH_FUNC == 5
  return 9 * pow(x, 10) - 100 * pow(x, 9);
#else
  return sin(x)/sin(1) - x;
#endif
}

#endif
