#ifndef CARTDG_CONNECTION_HPP_
#define CARTDG_CONNECTION_HPP_

#include "Deformed_element.hpp"

namespace cartdg
{

class Cartesian_face_connection
{
  public:
  virtual int i_dim() = 0;
  virtual double* face(int i_side) = 0;
};

class Deformed_face_connection
{
  public:
  virtual int i_dim(int i_side) = 0;
  virtual bool face_sign(int i_side) = 0;
  virtual double* face(int i_side) = 0;
  virtual double* jacobian() = 0;
};

class Element_connection
{
  public:
  virtual Element& element(int i_side) = 0;
};

class Cartesian_element_connection : public Cartesian_face_connection, public Element_connection
{
  std::array<Element*, 2> elems;

  public:
  Cartesian_element_connection(std::array<Element*, 2> elements, int i_dim_arg);
  virtual int i_dim();
  virtual double* face(int i_side);
  virtual Element& element(int i_side);
};

class Deformed_element_connection : public Deformed_face_connection, public Element_connection
{
  std::array<Deformed_element*, 2> elems;

  public:
  Deformed_element_connection(std::array<Deformed_element*, 2> elements, std::array<int, 2> i_dim_arg, std::array<bool, 2> face_sign_arg);
  virtual int i_dim(int i_side);
  virtual bool face_sign(int i_side);
  virtual double* face(int i_side);
  virtual double* jacobian();
  virtual Deformed_element& element(int i_side);
};

}
#endif
