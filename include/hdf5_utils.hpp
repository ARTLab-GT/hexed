#ifndef HEXED_HDF5_UTILS_HPP_
#define HEXED_HDF5_UTILS_HPP_

#include <H5Cpp.h>
#include "utils.hpp"

//! \brief Convenience wrappers for some common and verbose HDF5 operations.
namespace hexed::hdf5_utils {

//! \brief Fetches the HDF5 `DataType` object that describes the type `T`.
template <typename T>
const H5::PredType& type() {
  static_assert(always_false<T>(), "Not implemented for this type.");
}

template<> const H5::PredType& type<double>();
template<> const H5::PredType& type<Int>();
template<> const H5::PredType& type<int>();

//! \brief Adds a scalar attribute `attr_name` with value `value` to `obj`.
template <typename T>
void add_attr(H5::H5Object& obj, std::string attr_name, T value) {
  hsize_t attr_dim = 1;
  H5::DataSpace dspace(1, &attr_dim);
  auto attr = obj.createAttribute(attr_name, type<T>(), dspace);
  attr.write(type<T>(), &value);
}

//! \brief Returns the value of attribute `attr_name` from `obj`.
//! \details `attr_name` must have been previously created with C++ type `T`.
template <typename T>
T get_attr(H5::H5Object& obj, std::string attr_name) {
  auto attr = obj.openAttribute(attr_name.c_str());
  T value;
  attr.read(type<T>(), &value);
  return value;
}

//! \brief Returns the value of attribute `attr_name` from HDF5 file named `file_name`.
//! \details `attr_name` must have been previously created with C++ type `T`.
template <typename T>
T get_attr(std::string file_name, std::string attr_name) {
  H5::H5File file(file_name, H5F_ACC_RDONLY);
  return get_attr<T>(file, attr_name);
}

//! \brief Writes an entry to a 2D HDF5 DataSet.
template <typename T>
void write(H5::DataSet& dset, Int row, Int col, T value) {
  hsize_t dims [2] {1, 1};
  H5::DataSpace mspace (2, dims, nullptr);
  hsize_t offset [2] {hsize_t(row), hsize_t(col)};
  hsize_t stride [2] {1, 1};
  hsize_t block [2] {1, 1};
  auto dspace = dset.getSpace();
  dspace.selectHyperslab(H5S_SELECT_SET, dims, offset, stride, block);
  dset.write(&value, dset.getDataType(), mspace, dspace);
}

//! \brief Reads an entry from a 2D HDF5 DataSet.
template <typename T>
T read(H5::DataSet& dset, Int row, Int col) {
  hsize_t dims [2] {1, 1};
  H5::DataSpace mspace (2, dims, nullptr);
  hsize_t offset [2] {hsize_t(row), hsize_t(col)};
  hsize_t stride [2] {1, 1};
  hsize_t block [2] {1, 1};
  auto dspace = dset.getSpace();
  dspace.selectHyperslab(H5S_SELECT_SET, dims, offset, stride, block);
  T data;
  dset.read(&data, dset.getDataType(), mspace, dspace);
  return data;
}

}
#endif
