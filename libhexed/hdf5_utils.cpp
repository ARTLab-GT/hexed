#include <hexed/hdf5_utils.hpp>

namespace hexed::hdf5_utils {

template<> const H5::PredType& type<double>() {return H5::PredType::NATIVE_DOUBLE;}
template<> const H5::PredType& type<Int>() {return H5::PredType::NATIVE_LLONG;}
template<> const H5::PredType& type<int>() {return H5::PredType::NATIVE_INT;}

}
