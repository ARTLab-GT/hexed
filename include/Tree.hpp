#ifndef HEXED_TREE_HPP_
#define HEXED_TREE_HPP_

#include "math.hpp"

namespace hexed
{

class Tree
{
  Mat<> orig;
  double root_sz;
  int ref_level;
  Eigen::VectorXi coords;

  public:
  Tree(int n_dim, double root_size, Mat<> origin = Mat<>::Zero(3));

  //! \name basic instance information
  //!\{
  const int n_dim;
  Mat<> origin() const;
  int refinement_level() const;
  Eigen::VectorXi coordinates() const;
  double nominal_size() const;
  Mat<> nominal_position() const;
  //!\}

  //! \name parent/child status
  //!\{
  Tree* parent();
  std::vector<Tree*> children();
  inline bool is_root() const;
  inline bool is_leaf() const;
  //!\}

  //! \name modifiers
  //!\{
  void refine();
  void unrefine();
  //!\}

  //! \name traversing functions
  //!\{
  Tree* find_leaf(int ref_level, Eigen::VectorXi coords, Eigen::VectorXi bias = Eigen::VectorXi::Zero(3));
  Tree* find_leaf(Mat<> position, Eigen::VectorXi bias = Eigen::VectorXi::Zero(3));
  Tree* find_neighbor(Eigen::MatrixXi direction);
  //!\}
};

}
#endif
