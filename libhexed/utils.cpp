#include <hexed/utils.hpp>

namespace hexed {

Mat<> resize(const Mat<>& vec, Int size) {
  Mat<> resized = Mat<>::Zero(size);
  auto seq = Eigen::seqN(0, std::min<Int>(vec.size(), size));
  resized(seq) = vec(seq);
  return resized;
}

std::string file_extension(std::string file_name) {
  unsigned extension_start = file_name.find_last_of(".");
  HEXED_ASSERT(extension_start != std::string::npos, "`file_name` has no extension");
  std::string case_sensitive = file_name.substr(extension_start + 1, std::string::npos);
  std::string ext = case_sensitive;
  for (char& c : ext) c = tolower(c);
  return ext;
}

std::string to_string(int i) {return std::to_string(i);}
std::string to_string(Int i) {return std::to_string(i);}
std::string to_string(double d) {return format_str(100, "%+.6e", d);}
std::string to_string(std::string s) {return s;}
std::string to_string(bool b) {return b ? "true" : "false";}

std::string to_string(Mat<dyn, dyn> mat) {
  if (mat.cols() <= 1) {
    return "Mat<>{" + to_string(mat.data(), mat.size()) + "}";
  }
  std::string s = "Mat<>{\n";
  for (Int row = 0; row < (Int)mat.rows(); ++row) {
    s += "  ";
    for (Int col = 0; col < (Int)mat.cols(); ++col) {
      s += to_string(mat(row, col)) + ", ";
    }
    s.pop_back();
    s += "\n";
  }
  s += "}";
  return s;
}

}
