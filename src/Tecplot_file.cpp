#include <TECIO.h>
#include <MASTER.h>

#include <Tecplot_file.hpp>
#include <math.hpp>

namespace cartdg
{

Tecplot_file::Tecplot_file(std::string file_name, int n_dim, int n_var, int row_size, double time)
: n_dim{n_dim}, n_var{n_var}, row_size{row_size}, n_qpoint{custom_math::pow(row_size, n_dim)},
  time{time}, strand_id{1}, i_zone{0}
{
  INTEGER4 Debug = 0;
  INTEGER4 VIsDouble = 1;
  INTEGER4 FileType = 0;
  INTEGER4 fileFormat = 1; // 0 == PLT, 1 == SZPLT
  /*
  * Open the file and write the tecplot datafile
  * header information
  */

  std::string var_names = "";
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    var_names += " x" + std::to_string(i_dim);
  }
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    var_names += " q" + std::to_string(i_var);
  }
  var_names.erase(0, 1);
  TECINI142((char*)"flow solution", var_names.c_str(), file_name.c_str(), ".", &fileFormat, &FileType, &Debug, &VIsDouble);
}

Tecplot_file::~Tecplot_file()
{
  TECEND142();
}

void Tecplot_file::write_block(double* pos, double* vars)
{
  INTEGER4 ICellMax = 0;
  INTEGER4 JCellMax = 0;
  INTEGER4 KCellMax = 0;
  INTEGER4 DIsDouble = 1;
  INTEGER4 unused = 0; // ParentZone is no longer used
  INTEGER4 IsBlock = 1; /* Block */
  INTEGER4 NFConns = 0;
  INTEGER4 FNMode = 0;
  INTEGER4 TotalNumFaceNodes = 1;
  INTEGER4 TotalNumBndryFaces = 1;
  INTEGER4 TotalNumBndryConnections = 1;
  INTEGER4 ShrConn = 0;
  /*Ordered Zone Parameters*/
  INTEGER4 IMax = row_size;
  INTEGER4 JMax = (n_dim >= 2) ? row_size : 1;
  INTEGER4 KMax = (n_dim >= 3) ? row_size : 1;

  /* Ordered Zone */
  INTEGER4 ZoneType = 0;
  TECZNE142(("zone " + std::to_string(i_zone)).c_str(),
  &ZoneType,
  &IMax,
  &JMax,
  &KMax,
  &ICellMax,
  &JCellMax,
  &KCellMax,
  &time,
  &strand_id,
  &unused,
  &IsBlock,
  &NFConns,
  &FNMode,
  &TotalNumFaceNodes,
  &TotalNumBndryFaces,
  &TotalNumBndryConnections,
  NULL,
  NULL,
  NULL,
  &ShrConn);
  ++strand_id;
  ++i_zone;

  int size {n_dim*n_qpoint};
  TECDAT142(&size, pos, &DIsDouble);
  size = n_var*n_qpoint;
  TECDAT142(&size, vars, &DIsDouble);
}

void Tecplot_file::write_line_segment(double* pos, double* vars)
{
  INTEGER4 ZoneType {1};
  INTEGER4 NumPoints {row_size};
  INTEGER4 NumElements {row_size - 1};
  INTEGER4 NumFaces {0};
  INTEGER4 ICellMax {0};
  INTEGER4 JCellMax {0};
  INTEGER4 KCellMax {0};
  INTEGER4 ParentZone {0};
  INTEGER4 IsBlock {1};
  INTEGER4 NumFaceConnections {0};
  INTEGER4 FaceNeighborMode {0};
  INTEGER4 TotalNumFaceNodes {0};
  INTEGER4 NumConnectedBoundaryFaces {0};
  INTEGER4 TotalNumBoundaryConnections {0};
  INTEGER4* PassiveVarList {nullptr};
  INTEGER4* ValueLocation {nullptr};
  INTEGER4* SharVarFromZone {nullptr};
  INTEGER4 ShareConnectivityFromZone {0};
  TECZNE142(("zone " + std::to_string(i_zone)).c_str(),
            &ZoneType, &NumPoints, &NumElements, &NumFaces, &ICellMax, &JCellMax,
            &KCellMax, &time, &strand_id, &ParentZone, &IsBlock, &NumFaceConnections,
            &FaceNeighborMode, &TotalNumFaceNodes, &NumConnectedBoundaryFaces,
            &TotalNumBoundaryConnections, PassiveVarList, ValueLocation, SharVarFromZone,
            &ShareConnectivityFromZone);
  ++strand_id;
  ++i_zone;

  INTEGER4 IsDouble {1};
  INTEGER4 size {n_dim*row_size};
  TECDAT142(&size, pos, &IsDouble);
  size = n_var*row_size;
  TECDAT142(&size, vars, &IsDouble);

  INTEGER4 node_inds [2*(row_size - 1)];
  for (int i_elem = 0; i_elem < row_size - 1; ++i_elem)
  {
    node_inds[2*i_elem] = i_elem + 1;
    node_inds[2*i_elem + 1] = i_elem + 2;
  }
  TECNOD142(node_inds);
}

}
