#ifndef N_EXTREMA_HPP_
#define N_EXTREMA_HPP_

namespace cartdg
{
template<int row_size>
int n_extrema(double* data)
{
  int n = 0;
  bool positive_slope = (data[1] > data[0]);
  bool negative_slope = (data[1] < data[0]);
  for (int i = 1; i < row_size - 1; ++i)
  {
    bool new_positive_slope = (data[i + 1] > data[i]);
    bool new_negative_slope = (data[i + 1] < data[i]);
    bool is_max = positive_slope && new_negative_slope;
    bool is_min = negative_slope && new_positive_slope;
    if (is_max || is_min) ++n;
    positive_slope = new_positive_slope;
    negative_slope = new_negative_slope;
  }
  return n;
}

template<>
int n_extrema<0>(double* data) {return 0;}

template<>
int n_extrema<1>(double* data) {return 0;}

template<>
int n_extrema<2>(double* data) {return 0;}

}

#endif