#ifndef HEXED_ENTITY_HPP_
#define HEXED_ENTITY_HPP_

#include <memory>
#include <math.hpp>
#include <Mutual_ptr.hpp>
#include <Array.hpp>
#include <Basis.hpp>

namespace hexed
{

class Entity
{
  public:
  enum type {cartesian, interp, free};
  const int n_dim;
  const int my_dim;
  const std::shared_ptr<Basis> basis;
  Entity(int n_dim, int my_dim, const std::shared_ptr<Basis>);
  virtual type get_type() const = 0;
  virtual const Array<double> points() const = 0;
};

class Cartesian : public Entity
{
  Array<int> _nom_pos;
  public:
  Cartesian(int n_dim, int my_dim, const std::shared_ptr<Basis>, Array<int> nom_pos, double nom_sz);
  Array<int> nominal_position();
  inline type get_type() const {return cartesian;}
  const Array<double> points() const;
};

}
#endif
