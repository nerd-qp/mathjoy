/*
 * Implements functionality of the utils::function
 * Given its growing need of high performance
 */

#include "utils/utils.hpp"

namespace utils {

  typedef std::function<double (double, double)> func;

  const function &function::operator= (const func &f) {
    _self = f;
    return *this;
  }

  const function &function::operator= (const function &f) {
    this->_self = f.getFunction();
    return *this;
  }

  const function &function::operator= (func &&f) {
    _self = std::move(f);
    return *this;
  }

  const function &function::operator= (function &&f) {
    this->_self = std::move(f.getFunction());
    return *this;
  }

  // function(const func && f) {
  //   _self = f;
  // }

  double function::operator() (double x, double y) const {
    return _self(x, y);
  }

  function function::operator+ (const function &f) const {
    auto f1 = this->_self;
    auto f2 = f.getFunction();
    function res ([=] (double x, double y) {
                    return f1(x, y) + f2(x, y);
                  });
    return res;
  }

  function function::operator+ (double scalar) const {
    auto f1 = this->_self;
    function res([=] (double x, double y) -> double {
                   return f1(x, y) + scalar;
                 });
    return res;
  }

  const function &function::operator+= (const function &f) {
    *this = *this + f;
    return *this;
  }

  function function::operator* (const function &f) const {
    auto f1 = this->_self;
    auto f2 = f.getFunction();
    function res ([=] (double x, double y) {
                    return f1(x, y) * f2(x, y);
                  });
    return res;
  }

  function function::operator* (double scalar) const {
    auto f1 = this->_self;

    function res ([=] (double x, double y) -> double {
                    return f1(x, y) * scalar;
                  });
    return res;
  }

  function function::operator/ (const function &f) const {
    auto f1 = this->_self;
    auto f2 = f.getFunction();
    function res ([=] (double x, double y) {
                    return f1(x, y) / f2(x, y);
                  });
    return res;
  }

  void function::setFunction(const func &f) {
    this->_self = f;
  }

  void function::setFunction(func &&f) {
    this->_self = std::move(f);
  }

  const func &function::getFunction() const {
    return this->_self;
  }
}
