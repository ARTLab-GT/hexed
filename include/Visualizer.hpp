#ifndef HEXED_VISUALIZER_HPP_
#define HEXED_VISUALIZER_HPP_

namespace hexed
{

class Visualizer
{
  public:
  /*! \brief writes a structured block of data
   * \param row_size size of each row of sample points
   * \param pos layout: [i_dim][i_point]
   * \param vars layout: [i_var][i_point]
   */
  virtual void write_block(int row_size, double* pos, double* vars) = 0;
};

}
#endif
