#include "TECIO.h"
#include "MASTER.h"
#include <iostream>

#include "Grid.hpp"

void Grid::visualize(std::string file_name)
{
  INTEGER4 Debug = 1;
  INTEGER4 VIsDouble = 1;
  INTEGER4 FileType = 0;
  INTEGER4 fileFormat = 1; // 0 == PLT, 1 == SZPLT
  INTEGER4 I = 0; /* Used to track return codes */
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
  I = TECINI142((char*)"flow solution", var_names.c_str(), file_name.c_str(), ".",
                &fileFormat, &FileType, &Debug, &VIsDouble);

  INTEGER4 ICellMax = 0;
  INTEGER4 JCellMax = 0;
  INTEGER4 KCellMax = 0;
  INTEGER4 DIsDouble = 0;
  double SolTime = time;
  INTEGER4 StrandID = 0; /* StaticZone */
  INTEGER4 unused = 0; // ParentZone is no longer used
  INTEGER4 IsBlock = 1; /* Block */
  INTEGER4 NFConns = 0;
  INTEGER4 FNMode = 0;
  INTEGER4 TotalNumFaceNodes = 1;
  INTEGER4 TotalNumBndryFaces = 1;
  INTEGER4 TotalNumBndryConnections = 1;
  INTEGER4 ShrConn = 0;
  /*Ordered Zone Parameters*/
  INTEGER4 IMax = basis.rank;
  INTEGER4 JMax = (n_dim >= 2) ? basis.rank : 1;
  INTEGER4 KMax = (n_dim >= 3) ? basis.rank : 1;

  int i_elem = 0;
  /* Ordered Zone */
  INTEGER4 ZoneType = 0;
  I = TECZNE142((char*)"element 0",
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

  INTEGER4 III = IMax * JMax * KMax;
  std::vector<double> pos = get_pos(0);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    I = TECDAT142(&III, pos.data() + i_dim*n_qpoint, &DIsDouble);
  }
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    I = TECDAT142(&III, state_r.data() + i_elem*n_dof + i_var*n_qpoint, &DIsDouble);
  }

  I = TECEND142();
}
