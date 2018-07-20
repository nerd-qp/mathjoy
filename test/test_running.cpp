#include <gtest/gtest.h>
//#include <gmock/gmock.h>

// matrix and vector from boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "matrix_solver/gaussian_elimination.hpp"
#include "utils/matrix_io.hpp"


namespace {

  // // The fixture for testing class Foo.
  // class FooTest : public ::testing::Test {
  // protected:
  //   // You can remove any or all of the following functions if its body
  //   // is empty.

  //   FooTest() {
  //     // You can do set-up work for each test here.
  //   }

  //   ~FooTest() override {
  //     // You can do clean-up work that doesn't throw exceptions here.
  //   }

  //   // If the constructor and destructor are not enough for setting up
  //   // and cleaning up each test, you can define the following methods:

  //   void SetUp() override {
  //     // Code here will be called immediately after the constructor (right
  //     // before each test).
  //   }

  //   void TearDown() override {
  //     // Code here will be called immediately after each test (right
  //     // before the destructor).
  //   }

  //   // Objects declared here can be used by all tests in the test case for Foo.
  // };

  // // Tests that the Foo::Bar() method does Abc.
  // TEST_F(FooTest, MethodBarDoesAbc) {
  //   const std::string input_filepath = "this/package/testdata/myinputfile.dat";
  //   const std::string output_filepath = "this/package/testdata/myoutputfile.dat";
  //   Foo f;
  //   EXPECT_EQ(f.Bar(input_filepath, output_filepath), 0);
  // }

  // // Tests that Foo does Xyz.
  // TEST_F(FooTest, DoesXyz) {
  //   // Exercises the Xyz feature of Foo.
  // }

  using namespace boost::numeric::ublas;

  TEST(GaussianElimination, IdentityMatrixTest) {
    //using namespace boost::assign;
    matrix<double> A(identity_matrix<double>(3));
    //matrix<double> A(3, 3);
    vector<double> b(3);
    b(0) = b(1) = b(2) = 1;
    auto b_backup = b;
    auto A_backup = A;
    mathjoy::Gaussian_elimination(A, b);

    for ( auto i = 0ul; i < b_backup.size(); ++ i )
      ASSERT_DOUBLE_EQ(b(i), b_backup(i));

    for ( auto i = 0ul; i < A_backup.size1(); ++ i ) {
      for ( auto j = 0ul; j < A_backup.size2(); ++ j ) {
        ASSERT_DOUBLE_EQ(A(i, j), A_backup(i, j));
      }
    }

    //DoubleEq(1.0);
    //ASSERT_DOUBLE_EQ(1, 1);
    // ASSERT_THAT(b, )
    //ASSERT_EQ()

  }

  TEST(GaussianElimination, SolveMatrixOne) {
    auto A = utils::makeMatrix(3, 3,
                               std::vector<double>
                               {1, 0.5, 1.5,
                                2, 2, 6,
                                -1, -0.5, -0.5});

    vector<double, std::vector<double>> b({3.0, 10.0, -2.0});

    mathjoy::Gaussian_elimination(A, b);

    vector<double, std::vector<double>> expected_b({3.0, 4.0, 1.0});
    auto expected_A = utils::makeMatrix(3, 3,
                                        std::vector<double>
                                        {1, 0.5, 1.5,
                                         0, 1, 3,
                                         0, 0, 1});

    for ( auto i = 0ul; i < expected_b.size(); ++ i )
      ASSERT_DOUBLE_EQ(b(i), expected_b(i));

    for ( auto i = 0ul; i < expected_A.size1(); ++ i ) {
      for ( auto j = 0ul; j < expected_A.size2(); ++ j ) {
        ASSERT_DOUBLE_EQ(A(i, j), expected_A(i, j));
      }
    }

  }

}  // namespace

// just link gtest_main, magic
// int main(int argc, char **argv) {
//   ::testing::InitGoogleTest(&argc, argv);
//   return RUN_ALL_TESTS();
// }
