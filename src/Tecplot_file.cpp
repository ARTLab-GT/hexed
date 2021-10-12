#include <TECIO.h>
#include <MASTER.h>

#include <Tecplot_file.hpp>
#include <math.hpp>

namespace cartdg
{

int Tecplot_file::n_instances {0};

Tecplot_file::Tecplot_file(std::string file_name, int n_dim, int n_var, int row_size, double time)
: n_dim{n_dim}, n_var{n_var}, row_size{row_size}, n_qpoint{custom_math::pow(row_size, n_dim)},
  time{time}, strand_id{1}, i_zone{0}
{
  if (n_instances > 0) throw std::runtime_error("Attempt to create multiple `Tecplot_file`s at once, which is illegal.");
  ++n_instances;

  INTEGER4 Debug = 0;
  INTEGER4 VIsDouble = 1;
  INTEGER4 FileType = 0;
  INTEGER4 fileFormat = 1; // 0 indicates .plt and 1 indicates .szplt
  std::string var_names = "";
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    var_names += " position" + std::to_string(i_dim);
  }
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    var_names += " state" + std::to_string(i_var);
  }
  var_names.erase(0, 1);
  TECINI142((char*)"flow solution", var_names.c_str(), file_name.c_str(), ".", &fileFormat, &FileType, &Debug, &VIsDouble); // opens a new Tecplot file
}

Tecplot_file::~Tecplot_file()
{
  TECEND142(); // closes Tecplot file
  --n_instances;
}

void Tecplot_file::write_line_segment(double* pos, double* vars)
{
  INTEGER4 ZoneType {1}; // 1 indicates unstructured ("finite element", in Tecplot parlance) line segment
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
            &ShareConnectivityFromZone); // initialize new zone
  ++strand_id;
  ++i_zone;

  // write data
  INTEGER4 IsDouble {1};
  INTEGER4 size {n_dim*row_size};
  TECDAT142(&size, pos, &IsDouble);
  size = n_var*row_size;
  TECDAT142(&size, vars, &IsDouble);

  // define connections between nodes
  INTEGER4 node_inds [2*(row_size - 1)];
  for (int i_elem = 0; i_elem < row_size - 1; ++i_elem)
  {
    node_inds[2*i_elem] = i_elem + 1;
    node_inds[2*i_elem + 1] = i_elem + 2;
  }
  TECNOD142(node_inds);
}

Tecplot_file::Zone::Zone(Tecplot_file& file) : file{file} {}


void Tecplot_file::Zone::write(double* pos, double* vars)
{
  INTEGER4 IsDouble {1};
  INTEGER4 size {file.n_dim*n_nodes};
  TECDAT142(&size, pos, &IsDouble);
  size = file.n_var*n_nodes;
  TECDAT142(&size, vars, &IsDouble);
}

Tecplot_file::Structured_block::Structured_block(Tecplot_file& file)
: Zone{file}
{
  n_nodes = file.n_qpoint;

  INTEGER4 ICellMax = 0;
  INTEGER4 JCellMax = 0;
  INTEGER4 KCellMax = 0;
  INTEGER4 unused = 0; // ParentZone is no longer used
  INTEGER4 IsBlock = 1; // always set to 1
  INTEGER4 NFConns = 0;
  INTEGER4 FNMode = 0;
  INTEGER4 TotalNumFaceNodes = 1;
  INTEGER4 TotalNumBndryFaces = 1;
  INTEGER4 TotalNumBndryConnections = 1;
  INTEGER4 ShrConn = 0;
  INTEGER4 IMax = file.row_size;
  INTEGER4 JMax = (file.n_dim >= 2) ? file.row_size : 1;
  INTEGER4 KMax = (file.n_dim >= 3) ? file.row_size : 1;

  INTEGER4 ZoneType = 0; // 0 indicates ordered
  TECZNE142(("zone " + std::to_string(file.i_zone)).c_str(),
            &ZoneType, &IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax, &file.time, &file.strand_id,
            &unused, &IsBlock, &NFConns, &FNMode, &TotalNumFaceNodes, &TotalNumBndryFaces,
            &TotalNumBndryConnections, NULL, NULL, NULL, &ShrConn); // initialize new zone
  ++file.strand_id; // next zone will be on next time strand
  ++file.i_zone; // next zone will be named with next number
}

Tecplot_file::Line_segments::Line_segments(Tecplot_file& file, int n_segs)
: Zone{file}, n_segs{n_segs}, i_seg{0}
{
  n_nodes = file.row_size*n_segs;
  pos_storage.resize(file.n_dim*n_nodes);
  var_storage.resize(file.n_var*n_nodes);

  INTEGER4 ZoneType {1}; // 1 indicates unstructured ("finite element", in Tecplot parlance) line segment
  INTEGER4 NumPoints {n_nodes};
  INTEGER4 NumElements {(file.row_size - 1)*n_segs};
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
  TECZNE142(("zone " + std::to_string(file.i_zone)).c_str(),
            &ZoneType, &NumPoints, &NumElements, &NumFaces, &ICellMax, &JCellMax,
            &KCellMax, &file.time, &file.strand_id, &ParentZone, &IsBlock, &NumFaceConnections,
            &FaceNeighborMode, &TotalNumFaceNodes, &NumConnectedBoundaryFaces,
            &TotalNumBoundaryConnections, PassiveVarList, ValueLocation, SharVarFromZone,
            &ShareConnectivityFromZone); // initialize new zone
  ++file.strand_id;
  ++file.i_zone;
}

void Tecplot_file::Line_segments::write(double* pos, double* vars)
{
  for (int i_node = 0; i_node < file.row_size; ++i_node)
  {
    for (int i_dim = 0; i_dim < file.n_dim; ++i_dim)
    {
      pos_storage[(n_segs*i_dim + i_seg)*file.row_size + i_node] = pos[i_dim*file.row_size + i_node];
    }
    for (int i_var = 0; i_var < file.n_var; ++i_var)
    {
      var_storage[(n_segs*i_var + i_seg)*file.row_size + i_node] = vars[i_var*file.row_size + i_node];
    }
  }
  ++i_seg;
}

Tecplot_file::Line_segments::~Line_segments()
{
  Zone::write(pos_storage.data(), var_storage.data());
  std::vector<INTEGER4> inds;
  for (int i_seg = 0; i_seg < n_segs; ++i_seg)
  {
    for (int i_elem = 0; i_elem < file.row_size - 1; ++i_elem)
    {
      inds.push_back(file.row_size*i_seg + i_elem + 1);
      inds.push_back(file.row_size*i_seg + i_elem + 2);
    }
  }
  TECNOD142(inds.data());
}

}
