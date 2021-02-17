
/*
This file was generated automatically by script/auto_generate.py.
Do not attempt to modify it directly. Instead, modify and rerun script/auto_generate.py
to make the required changes.
*/


#include <Solution.hpp>
#include <kernels/local/cpg_euler_matrix.hpp>
#include <kernels/neighbor/cpg_euler_copy.hpp>
#include <kernels/neighbor/cpg_euler_fbc.hpp>
#include <kernels/max_char_speed/cpg_euler_max.hpp>

namespace cartdg
{

Local_kernel local_kernels [3][5] {
        &(cpg_euler_matrix<3, 1, 1>),
                &(cpg_euler_matrix<3, 2, 2>),
                &(cpg_euler_matrix<3, 3, 3>),
                &(cpg_euler_matrix<3, 4, 4>),
                &(cpg_euler_matrix<3, 5, 5>),
                &(cpg_euler_matrix<4, 1, 1>),
                &(cpg_euler_matrix<4, 4, 2>),
                &(cpg_euler_matrix<4, 9, 3>),
                &(cpg_euler_matrix<4, 16, 4>),
                &(cpg_euler_matrix<4, 25, 5>),
                &(cpg_euler_matrix<5, 1, 1>),
                &(cpg_euler_matrix<5, 8, 2>),
                &(cpg_euler_matrix<5, 27, 3>),
                &(cpg_euler_matrix<5, 64, 4>),
                &(cpg_euler_matrix<5, 125, 5>),
        };

    Local_kernel Solution::get_local_kernel()
    {
      if (n_dim <= 3) return local_kernels[n_dim - 1][basis.rank - 1];
      else throw "Kernel not available.";
    }
    
Neighbor_kernel neighbor_kernels [3][5] {
        &(cpg_euler_copy<3, 1, 1>),
                &(cpg_euler_copy<3, 2, 2>),
                &(cpg_euler_copy<3, 3, 3>),
                &(cpg_euler_copy<3, 4, 4>),
                &(cpg_euler_copy<3, 5, 5>),
                &(cpg_euler_copy<4, 1, 1>),
                &(cpg_euler_copy<4, 4, 2>),
                &(cpg_euler_copy<4, 9, 3>),
                &(cpg_euler_copy<4, 16, 4>),
                &(cpg_euler_copy<4, 25, 5>),
                &(cpg_euler_copy<5, 1, 1>),
                &(cpg_euler_copy<5, 8, 2>),
                &(cpg_euler_copy<5, 27, 3>),
                &(cpg_euler_copy<5, 64, 4>),
                &(cpg_euler_copy<5, 125, 5>),
        };

    Neighbor_kernel Solution::get_neighbor_kernel()
    {
      if (n_dim <= 3) return neighbor_kernels[n_dim - 1][basis.rank - 1];
      else throw "Kernel not available.";
    }
    
Fbc_kernel fbc_kernels [3][5] {
        &(cpg_euler_fbc<3, 1, 1>),
                &(cpg_euler_fbc<3, 2, 2>),
                &(cpg_euler_fbc<3, 3, 3>),
                &(cpg_euler_fbc<3, 4, 4>),
                &(cpg_euler_fbc<3, 5, 5>),
                &(cpg_euler_fbc<4, 1, 1>),
                &(cpg_euler_fbc<4, 4, 2>),
                &(cpg_euler_fbc<4, 9, 3>),
                &(cpg_euler_fbc<4, 16, 4>),
                &(cpg_euler_fbc<4, 25, 5>),
                &(cpg_euler_fbc<5, 1, 1>),
                &(cpg_euler_fbc<5, 8, 2>),
                &(cpg_euler_fbc<5, 27, 3>),
                &(cpg_euler_fbc<5, 64, 4>),
                &(cpg_euler_fbc<5, 125, 5>),
        };

    Fbc_kernel Solution::get_fbc_kernel()
    {
      if (n_dim <= 3) return fbc_kernels[n_dim - 1][basis.rank - 1];
      else throw "Kernel not available.";
    }
    
Max_char_speed_kernel max_char_speed_kernels [3][5] {
        &(cpg_euler_max<3, 1, 1>),
                &(cpg_euler_max<3, 2, 2>),
                &(cpg_euler_max<3, 3, 3>),
                &(cpg_euler_max<3, 4, 4>),
                &(cpg_euler_max<3, 5, 5>),
                &(cpg_euler_max<4, 1, 1>),
                &(cpg_euler_max<4, 4, 2>),
                &(cpg_euler_max<4, 9, 3>),
                &(cpg_euler_max<4, 16, 4>),
                &(cpg_euler_max<4, 25, 5>),
                &(cpg_euler_max<5, 1, 1>),
                &(cpg_euler_max<5, 8, 2>),
                &(cpg_euler_max<5, 27, 3>),
                &(cpg_euler_max<5, 64, 4>),
                &(cpg_euler_max<5, 125, 5>),
        };

    Max_char_speed_kernel Solution::get_max_char_speed_kernel()
    {
      if (n_dim <= 3) return max_char_speed_kernels[n_dim - 1][basis.rank - 1];
      else throw "Kernel not available.";
    }
    

}
