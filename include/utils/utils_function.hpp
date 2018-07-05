/*
 * std::function class wrapper to overload some operator
 */

/*
  Include libraries
 */

#include <functional>
#include <boost/numeric/ublas/vector.hpp>

namespace utils {

  // A very light proxy class of std::function
  // To support addition and multiplication semantics(operations)
  class function {
  private:
    typedef std::function<double (double, double)> func;
    func _self;

  public:
    // Maybe some copy construtor, move semantics
    // TODO Dig into this
    function() { };

    function(const func & f) {
      _self = f;
    }

    // template<typename Functor>
    // function(Functor functor) {
    //   _self = func(functor);
    // }

    function(const function &f) {
      //std::cout << "Copy!\n";
      this->_self = f.getFunction();
    }

    // function(function &&f) {
    //   std::cout << "Moved!\n";
    //   this->_self = std::move(f.getFunction());
    // }

    const function &operator= (const func &f);

    const function &operator= (const function &f);

    const function &operator= (func &&f);

    const function &operator= (function &&f);

    /*
     * operator(), functor capacity
     */
    double operator() (double x, double y) const;

    /*
     * operation between functions to return function
     * (temporary, please don't use const reference,
     *  it's just stupid)
     */
    function operator+ (const function &f) const;

    /*
     * Ability to add a scalar
     */
    function operator+ (double scalar) const;

    /*
     * Complete the operator+ above
     */
    friend function operator+ (double scalar, const function &f);

    const function &operator+= (const function &f);

    function operator* (const function &f) const;

    function operator* (double scalar) const;

    function operator/ (const function &f) const;

    void setFunction(const func &f);

    void setFunction(func &&f);

    const func &getFunction() const;
  };

  class function_1dim {
  private:
    typedef std::function<double (double)> func;
    func _self;

  public:
    function_1dim() { };

    function_1dim(const func &f) { _self = f; };

    function_1dim(const function_1dim &f) { _self = f.getFunction(); };

    function_1dim(func &&f) {
      _self = f;
    }

    function_1dim(function_1dim &&f) {
      _self = std::move(f.getFunction());
    }

    const function_1dim & operator= (const func &f) {
      _self = f;
    }

    const function_1dim & operator= (const function_1dim &f) {
      _self = f.getFunction();
    }

    double operator() (double x) const {
      return _self(x);
    }

    function operator* (const function_1dim &f) const {
      // make a copy of the two functions
      auto f1 = this->_self;
      auto f2 = f.getFunction();

      // NOTE Must avoid (*this) in lambda
      // It's dangerous, since that reference could be gone away
      // (= would capture this by reference for what??
      // does not make sense to capture by reference)
      function res ([=] (double x, double y) -> double {
                      // std::cout << x << ", " << y << '\n';
                      // std::cout << f1(x) * f2(y) << std::endl;
                      return f1(x) * f2(y);
                    });
      return res;
    }

    const func &getFunction() const {
      return _self;
    }
  };

  template<typename Vector>
  boost::numeric::ublas::vector<function>
  tensor_product_fun(const Vector &a,
                     const Vector &b) {
    boost::numeric::ublas::vector<function> res(a.size() * b.size());
    for ( size_t i = 0; i < a.size(); ++ i )
      for ( size_t j = 0; j < b.size(); ++ j )
        res(i * a.size() + j) = (a(i) * b(j));
    return res;
  }

  template<typename Vector>
  boost::numeric::ublas::vector<function>
  element_prod(const Vector &a,
               const Vector &b) {
    // a and b should have the same size
    if ( a.size() != b.size() )
      throw "Not the same size!\n";

    boost::numeric::ublas::vector<function> res(a.size());
    for ( size_t i = 0; i < a.size(); ++ i )
      res(i)= a(i) * b(i);

    return res;
  }
}
