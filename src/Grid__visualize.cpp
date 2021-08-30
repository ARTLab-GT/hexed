#include <iostream>

#include <TECIO.h>
#include <MASTER.h>

#include "Grid.hpp"

namespace cartdg
{

void Grid::visualize(std::string file_name)
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
  TECINI142((char*)"flow solution", var_names.c_str(), file_name.c_str(), ".",
            &fileFormat, &FileType, &Debug, &VIsDouble);

  INTEGER4 ICellMax = 0;
  INTEGER4 JCellMax = 0;
  INTEGER4 KCellMax = 0;
  INTEGER4 DIsDouble = 1;
  double SolTime = time;
  INTEGER4 StrandID = 1; // This must not be 0, as 0 indicats a StaticZone (apparently)
  INTEGER4 unused = 0; // ParentZone is no longer used
  INTEGER4 IsBlock = 1; /* Block */
  INTEGER4 NFConns = 0;
  INTEGER4 FNMode = 0;
  INTEGER4 TotalNumFaceNodes = 1;
  INTEGER4 TotalNumBndryFaces = 1;
  INTEGER4 TotalNumBndryConnections = 1;
  INTEGER4 ShrConn = 0;
  /*Ordered Zone Parameters*/
  INTEGER4 IMax = basis.row_size;
  INTEGER4 JMax = (n_dim >= 2) ? basis.row_size : 1;
  INTEGER4 KMax = (n_dim >= 3) ? basis.row_size : 1;

  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    /* Ordered Zone */
    INTEGER4 ZoneType = 0;
    TECZNE142(("element " + std::to_string(i_elem)).c_str(),
    &ZoneType,
    &IMax,
    &JMax,
    &KMax,
    &ICellMax,
    &JCellMax,
    &KCellMax,
    &SolTime,
    &StrandID,
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
    ++StrandID;

    INTEGER4 III = IMax * JMax * KMax;
    std::vector<double> pos = get_pos(i_elem);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      TECDAT142(&III, pos.data() + i_dim*n_qpoint, &DIsDouble);
    }
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      TECDAT142(&III, state_r() + i_elem*n_dof + i_var*n_qpoint, &DIsDouble);
    }
  }

  TECEND142();
}

}
