#include "catch.hpp"
#include "kernels/local/cpg_euler_matrix.hpp"
#include "kernels/local/cpg_euler_tensor.hpp"

#define TEST_CPG_EULER(kernel, cfl) \
  SECTION("1D") \
  { \
    const int n_elem = 5; \
    const int rank = 2; \
    double read [n_elem][3][rank]; \
    double write[n_elem][3][rank]; \
    Eigen::MatrixXd diff_mat = Eigen::MatrixXd::Identity(rank, rank); \
    double mass = 1.225; double veloc = 10; double pres = 1e5; \
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc; \
    for (int i_elem = 0; i_elem < n_elem; ++i_elem) \
    { \
      for (int i_qpoint = 0; i_qpoint < rank; ++i_qpoint) \
      { \
          read[i_elem][0][i_qpoint] = mmtm; \
          read[i_elem][1][i_qpoint] = mass; \
          read[i_elem][2][i_qpoint] = ener; \
      } \
    } \
    kernel<3, 2, 2>(&read[0][0][0], &write[0][0][0], n_elem, \
                              diff_mat, cfl); \
    for (int i_elem = 0; i_elem < n_elem; ++i_elem) \
    { \
      for (int i_qpoint = 0; i_qpoint < rank; ++i_qpoint) \
      { \
        CHECK((write[i_elem][0][i_qpoint] - read[i_elem][0][i_qpoint]) \
              == Approx(0.1*(mmtm*mmtm/mass + pres))); \
        CHECK((write[i_elem][1][i_qpoint] - read[i_elem][1][i_qpoint]) \
              == Approx(0.1*mmtm)); \
        CHECK((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint]) \
              == Approx(0.1*((ener + pres)*veloc))); \
      } \
    } \
  } \
  \
  SECTION("3D") \
  { \
    const int n_elem = 5; \
    const int rank = 3; \
    double read [n_elem][5][rank*rank*rank]; \
    double write[n_elem][5][rank*rank*rank]; \
    Eigen::MatrixXd diff_mat = Eigen::MatrixXd::Identity(rank, rank); \
    double mass = 1.225; \
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30; \
    double pres = 1e5; \
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2); \
    for (int i_elem = 0; i_elem < n_elem; ++i_elem) \
    { \
      for (int i_qpoint = 0; i_qpoint < rank*rank*rank; ++i_qpoint) \
      { \
          read[i_elem][0][i_qpoint] = mass*veloc0; \
          read[i_elem][1][i_qpoint] = mass*veloc1; \
          read[i_elem][2][i_qpoint] = mass*veloc2; \
          read[i_elem][3][i_qpoint] = mass; \
          read[i_elem][4][i_qpoint] = ener; \
      } \
    } \
    kernel<5, 27, 3>(&read[0][0][0], &write[0][0][0], n_elem, \
                               diff_mat, cfl); \
    for (int i_elem = 0; i_elem < n_elem; ++i_elem) \
    { \
      for (int i_qpoint = 0; i_qpoint < rank*rank*rank; ++i_qpoint) \
      { \
        CHECK((write[i_elem][0][i_qpoint] - read[i_elem][0][i_qpoint]) \
              == Approx(0.1*(mass*veloc0*(veloc0 + veloc1 + veloc2) + pres))); \
        CHECK((write[i_elem][1][i_qpoint] - read[i_elem][1][i_qpoint]) \
              == Approx(0.1*(mass*veloc1*(veloc0 + veloc1 + veloc2) + pres))); \
        CHECK((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint]) \
              == Approx(0.1*(mass*veloc2*(veloc0 + veloc1 + veloc2) + pres))); \
        CHECK((write[i_elem][3][i_qpoint] - read[i_elem][3][i_qpoint]) \
              == Approx(0.1*mass*(veloc0 + veloc1 + veloc2))); \
        CHECK((write[i_elem][4][i_qpoint] - read[i_elem][4][i_qpoint]) \
              == Approx(0.1*((ener + pres)*(veloc0 + veloc1 + veloc2)))); \
      } \
    } \
  } \

TEST_CASE("Local kernels")
{
  SECTION("CPG Euler flux, matrix form")
  {
    TEST_CPG_EULER(cpg_euler_matrix, 0.1)
  }

  /*
  SECTION("CPG Euler flux, tensor form")
  {
    TEST_CPG_EULER(cpg_euler_tensor, 0.05)
  }
  */
}
