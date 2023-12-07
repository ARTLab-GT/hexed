#ifndef HEXED_XDMF_HPP_
#define HEXED_XDMF_HPP_

#include "config.hpp"
#if HEXED_USE_XDMF
#include <XdmfAttribute.hpp>
#include <XdmfTopology.hpp>
#include <XdmfGeometry.hpp>
#include <Eigen/Dense>
#include "Output_data.hpp"
#include "Visualizer.hpp"

namespace hexed
{

//! \brief lightweight wrapper for [XDMF](https://www.xdmf.org/index.php/XDMF_Model_and_Format)
//! [API](https://www.xdmf.org/index.php/Xdmf3_C%2B%2B_API) for block-structured data
class Xdmf_wrapper : public Visualizer
{
  boost::shared_ptr<XdmfTopology> _topo;
  boost::shared_ptr<XdmfGeometry> _geom;
  std::vector<boost::shared_ptr<XdmfAttribute>> _attrs;
  const int _n_dim_geom;
  const int _n_dim_topo;
  const std::string _file_name;
  const double _time;
  const int _n_var;
  int _i_block = 0;
  Eigen::MatrixXi _node_inds;
  public:
  /*!
   * \param n_dim_geom Number of geometric dimensions.
   *   Distinct from the number of topological dimensions.
   *   For example, a 3D surface will have 3 geometric dimensions and 2 topological dimensions.
   * \param n_dim_topo defines the topology type/dimensionality
   * \param row_size row size of data sample blocks
   * \param file_name name of output file(s) without extension
   * \param data defines the number and names of variables to be visualized
   * \param time flow time
   */
  Xdmf_wrapper(int n_dim_geom, int n_dim_topo, std::string file_name, const Output_data& data, double time);
  void write_block(int row_size, double* pos, double* vars) override;
  ~Xdmf_wrapper(); //!< writes the data to the file(s)
};

}
#endif
#endif
