#include <catch2/catch.hpp>

#include <Kernel_settings.hpp>
#include <Gauss_lobatto.hpp>
#include <get_nonpen_convective.hpp>

TEST_CASE("non-penetration BC")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;
  cartdg::Gauss_lobatto basis (2);

  double veloc0 = 100.;
  double veloc1 = 50.;
  double veloc2 = 0.;
  double mass = 1.2;
  double pressure = 1e5;
  double kin_ener = mass*0.5*(veloc0*veloc0 + veloc1*veloc1);
  double ener = (pressure/(settings.cpg_heat_rat - 1.) + kin_ener);

  cartdg::Storage_params params {3, 5, 3, 2};
  cartdg::def_elem_vec elements;
  elements.emplace_back(new cartdg::Deformed_element {params});
  elements.emplace_back(new cartdg::Deformed_element {params});
  cartdg::def_elem_wall_vec walls {{elements[1].get(), 1, 0},
                                   {elements[0].get(), 0, 1},
                                   {elements[1].get(), 1, 1}};

  for (int i_elem : {0, 1})
  {
    double* read = elements[i_elem]->stage(0);
    double* write = elements[i_elem]->stage(1);
    double* jacobian = elements[i_elem]->jacobian();
    for (int i_qpoint = 0; i_qpoint < 8; ++i_qpoint)
    {
      read[0*8 + i_qpoint] = veloc0*mass;
      read[1*8 + i_qpoint] = veloc1*mass;
      read[2*8 + i_qpoint] = veloc2*mass;
      read[3*8 + i_qpoint] = mass;
      read[4*8 + i_qpoint] = ener;
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        write[i_var*8 + i_qpoint] = 0.;
      }

      for (int i_jac = 0; i_jac < 9; ++i_jac) jacobian[i_jac*8 + i_qpoint] = 0.;
      jacobian[(2*3 + 0)*8 + i_qpoint] = 1;
      jacobian[(0*3 + 1)*8 + i_qpoint] = 1;
      jacobian[(1*3 + 2)*8 + i_qpoint] = 2;
      jacobian[(0*3 + 2)*8 + i_qpoint] = 1;
    }
  }

  cartdg::get_nonpen_convective(3, 2)(walls, basis, settings);
  {
    double* write = elements[0]->stage(1);
    for (int i_qpoint = 0; i_qpoint < 8; ++i_qpoint)
    {
      REQUIRE(write[0*8 + i_qpoint] == Approx(0.));
      REQUIRE(write[1*8 + i_qpoint] == Approx(0.));
      REQUIRE(write[2*8 + i_qpoint] == Approx(0.));
      REQUIRE(write[3*8 + i_qpoint] == Approx(0.));
      REQUIRE(write[4*8 + i_qpoint] == Approx(0.));
    }
  }
  {
    double* write = elements[1]->stage(1);
    double veloc_normal = 2*veloc0 - 1*veloc1;
    double mult = settings.d_t_by_d_pos/1.;
    double mom_tang = mass*(veloc0*1 + veloc1*2);
    for (int i_qpoint : {0, 1, 4, 5})
    {
      double d_mom_tang = (write[0*8 + i_qpoint]*1 + write[1*8 + i_qpoint]*2);
      REQUIRE(d_mom_tang == Approx(-veloc_normal*mom_tang*mult));
      REQUIRE(write[2*8 + i_qpoint] == Approx(0.));
      REQUIRE(write[3*8 + i_qpoint] == Approx(-veloc_normal*mass*mult));
      REQUIRE(write[4*8 + i_qpoint] == Approx(-veloc_normal*(ener + pressure)*mult));
    }
    for (int i_qpoint : {2, 3, 6, 7})
    {
      double d_mom_tang = (write[0*8 + i_qpoint]*1 + write[1*8 + i_qpoint]*2);
      REQUIRE(d_mom_tang == Approx(veloc_normal*mom_tang*mult));
      REQUIRE(write[2*8 + i_qpoint] == Approx(0.));
      REQUIRE(write[3*8 + i_qpoint] == Approx(veloc_normal*mass*mult));
      REQUIRE(write[4*8 + i_qpoint] == Approx(veloc_normal*(ener + pressure)*mult));
    }
  }
}
