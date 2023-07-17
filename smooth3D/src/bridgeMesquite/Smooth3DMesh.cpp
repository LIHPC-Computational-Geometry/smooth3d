/*
 * Smooth3DMesh.cpp
 *
 *  Created on: 21 févr. 2013
 *      Author: weilljc
 */

#include "Smooth3DMesh.h"

#include <cstddef>
#include <vector>

#include <MeshImplData.hpp>
#include <MsqError.hpp>
#include <PatchData.hpp>

Smooth3DMesh::Smooth3DMesh() {
  // TODO Auto-generated constructor stub

}

Smooth3DMesh::~Smooth3DMesh() {
  // TODO Auto-generated destructor stub
}

void Smooth3DMesh::initSmooth3DMesh(int_type nb_cells, int_type nb_nodes,
    const int_type* nb_nodes_per_cell, const int_type* nodes_number, double* x,
    double* y, double* z) {

  Mesquite2::MsqError err;

  numCoords = 3;
  clear();
  myMesh->allocate_vertices(nb_nodes, err);
  for (int_type vertex_id = 0; vertex_id < nb_nodes; ++vertex_id) {
    myMesh->reset_vertex(vertex_id,
	Mesquite2::Vector3D(x[vertex_id], y[vertex_id], z[vertex_id]), false,
	err);
  }

  myMesh->allocate_elements(nb_cells, err);
  std::vector<std::size_t> vertices;
  Mesquite2::EntityTopology elem_type;

  const int_type * ptr_connectivity = nodes_number;
  for (std::size_t cell_id = 0; cell_id < nb_cells; ++cell_id) {
    vertices.resize(nb_nodes_per_cell[cell_id]);
    for (size_t j = 0; j < nb_nodes_per_cell[cell_id]; j++)
      vertices[j] = *ptr_connectivity++;

    switch (nb_nodes_per_cell[cell_id]) {
    case 4: /* Tetra*/
      elem_type = Mesquite2::TETRAHEDRON;
      break;
    case 5: /* pyramide */
      elem_type = Mesquite2::PYRAMID;
      break;
    case 6: /* prisme */
      elem_type = Mesquite2::PRISM;
      break;
    case 8: /* Hexa */
      elem_type = Mesquite2::HEXAHEDRON;
      break;
    default:
      assert(0 == 1);
      break;
    }
    myMesh->reset_element(cell_id, vertices, elem_type, err);

  }

  mark_skin_fixed(err, true);
}

void Smooth3DMesh::epilogue(double* x, double* y, double* z,
    const double* relax) {

  Mesquite2::MsqError err;
  /* N'epiloguons pas ! */
  Mesquite2::VertexIterator* iter = vertex_iterator(err);
  int vertex_id = 0;
  while (!iter->is_at_end()) {
    Mesquite2::Mesh::VertexHandle handle = iter->operator*();
    iter->operator++();
    Mesquite2::MsqVertex coordinates;
    vertices_get_coordinates(&handle, &coordinates, 1, err);
    double rel = relax[vertex_id];
    if (rel != 0.0) {
      double m_x0, m_y0, m_z0;
      m_x0 = coordinates.x();
      m_y0 = coordinates.y();
      m_z0 = coordinates.z();
      x[vertex_id] = (1 - rel) * x[vertex_id] + rel * (m_x0);
      y[vertex_id] = (1 - rel) * y[vertex_id] + rel * (m_y0);
      z[vertex_id] = (1 - rel) * z[vertex_id] + rel * (m_z0);
    }
    ++vertex_id;
  }

}

/* les noeuds à relax == 0.0 sont fixed */
void Smooth3DMesh::markNonRelaxedFixed(const double* relax) {
  Mesquite2::MsqError err;
  for (size_t i = 0; i < myMesh->max_vertex_index(); ++i) {
    if (relax[i] == 0.0)
      myMesh->fix_vertex(i, true, err);
  }
}
