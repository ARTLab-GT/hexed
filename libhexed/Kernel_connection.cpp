#include <hexed/Kernel_connection.hpp>

namespace hexed {

bool operator==(Connection_direction dir0, Connection_direction dir1) {
  return dir0.i_dim == dir1.i_dim && dir0.face_sign == dir1.face_sign && dir0.rotate == dir1.rotate;
}

std::string to_string(Connection_direction dir) {
  return format_str("Connection_direction{i_dim = %i,%i; face_sign = %i,%i; rotate = %i}",
                    dir.i_dim[0], dir.i_dim[1], int(dir.face_sign[0]), int(dir.face_sign[1]), dir.rotate);
}

std::ostream& operator<<(std::ostream& stream, Connection_direction dir) {return stream << to_string(dir);}

}
