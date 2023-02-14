#include "pde.hpp"

namespace hexed::pde
{

std::shared_ptr<Transport_model> inviscid()
{
  return std::shared_ptr<Transport_model>(new Const_transport(0.));
}

std::shared_ptr<Transport_model> air_const_dyn_visc()
{
  return std::shared_ptr<Transport_model>(new Const_transport(1.81206e-5));
}

std::shared_ptr<Transport_model> air_const_therm_cond()
{
  return std::shared_ptr<Transport_model>(new Const_transport(1.81206e-5*1.4*287.0528/.4/.71));
}

std::shared_ptr<Transport_model> air_sutherland_dyn_visc()
{
  return std::shared_ptr<Transport_model>(new Sutherland(1.716e-5, 273., 111.));
}

std::shared_ptr<Transport_model> air_sutherland_therm_cond()
{
  return std::shared_ptr<Transport_model>(new Sutherland(.0241, 273., 194.));
}

}
