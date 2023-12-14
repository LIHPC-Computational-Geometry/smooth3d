/*
 * GetMe2D.cpp
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <cmath>
#include <cassert>
#include <cstdlib>		// pour rand

#include "smooth3D/smooth.h"
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>

#include "math/Real3.h"
#include "math/MeanRatio.h"

/*
// Fonction sp�cialis�e pour calculer le meanRatio pour un triangle
static Smooth3D::Real computeMeanRatio(gmds::Node& n0, gmds::Node& n1,
    gmds::Node& n2) {
  Smooth3D::Real3 p[3];
  p[0] = Smooth3D::Real3(n0.X(), n0.Y(), n0.Z());
  p[1] = Smooth3D::Real3(n1.X(), n1.Y(), n1.Z());
  p[2] = Smooth3D::Real3(n2.X(), n2.Y(), n2.Z());

  Smooth3D::Real3 normal = math::vecMul(p[1] - p[0], p[2] - p[0]);
  Smooth3D::Real norm = normal.abs();
  if (norm == 0.0)
    return (0.0);

  normal = normal / norm;

  Smooth3D::Real mean_ratio;

  Smooth3D::mFvn2e(mean_ratio, p, normal);

  return (mean_ratio);
}

// Fonction sp�cialis�e pour calculer le meanRatio pour un quadrilat�re
static Smooth3D::Real computeMeanRatio(gmds::Node& n0, gmds::Node& n1,
    gmds::Node& n2, gmds::Node& n3) {
  static const int kTriQuads[4][3] = { { 0, 1, 3 }, { 1, 2, 0 }, { 2, 3, 1 }, {
      3, 0, 2 } };
  Smooth3D::Real3 d_con(-1.0, -1.0, -1.0);

  Smooth3D::Real3 p[4];
  p[0] = Smooth3D::Real3(n0.X(), n0.Y(), n0.Z());
  p[1] = Smooth3D::Real3(n1.X(), n1.Y(), n1.Z());
  p[2] = Smooth3D::Real3(n2.X(), n2.Y(), n2.Z());
  p[3] = Smooth3D::Real3(n3.X(), n3.Y(), n3.Z());

  Smooth3D::Real3 normal = math::vecMul(p[2] - p[0], p[3] - p[1]);

  Smooth3D::Real norm = normal.abs();

  normal = normal / norm;
  Smooth3D::Real somme = 0.0;

  for (int i = 0; i < 4; i++) {
    Smooth3D::Real3 q[3];
    Smooth3D::Real m;

    q[0] = p[kTriQuads[i][0]];
    q[1] = p[kTriQuads[i][1]];
    q[2] = p[kTriQuads[i][2]];

    Smooth3D::mFcn2i(m, q, normal, d_con);
    somme += m;

  }
  return (0.25 * somme);
}

static Smooth3D::Real computeMeanRatio(gmds::Face& face) {
  std::vector<gmds::Node> nodes = face.get<gmds::Node>();
  switch (nodes.size()) {
  case 3:
    // Triangle
    return computeMeanRatio(nodes[0], nodes[1], nodes[2]);
  case 4:
    // Quadrangle
    return computeMeanRatio(nodes[0], nodes[1], nodes[2], nodes[3]);
  default:
    break;
  }

  return (0.0);
}

class FourNodes {
  Smooth3D::Real3 m_a[4];
public:
  FourNodes() {
    m_a[0] = m_a[1] = m_a[2] = m_a[3] = Smooth3D::Real3::null();
  }

  Smooth3D::Real3 & operator[](int idx) {
    return m_a[idx];
  }
};



static FourNodes computeGetMeTriangle3D(gmds::Face& face, Smooth3D::Real theta,
    Smooth3D::Real lambda) {
  FourNodes compute_nodes, new_nodes;

  std::vector<gmds::Node> triangle_nodes = face.get<gmds::Node>();
  for (int i = 0; i < 3; i++) {
    gmds::Node current_node = triangle_nodes[i];
    compute_nodes[i] = Smooth3D::Real3(current_node.X(), current_node.Y(),
	current_node.Z());
    // We map the nodes in 2D, Z0 => (0.0)

  }
  Smooth3D::Real3 centroid = (compute_nodes[0] + compute_nodes[1]
      + compute_nodes[2]) * (1.0 / 3.0);
  Smooth3D::Real perimeter = (compute_nodes[1] - compute_nodes[0]).abs()
      + (compute_nodes[2] - compute_nodes[1]).abs()
      + (compute_nodes[2] - compute_nodes[0]).abs();

  Smooth3D::Real3 normal = math::vecMul(compute_nodes[1] - compute_nodes[0],
      compute_nodes[2] - compute_nodes[1]);
  normal = normal / normal.abs();

  Smooth3D::Real t = std::tan(theta);

  Smooth3D::Real3 z0, z1, z2;
  // premiere passe

  for (int i = 0; i < 3; i++) {
    Smooth3D::Real3 z0 = compute_nodes[i];
    Smooth3D::Real3 z1 = compute_nodes[(i + 1) % 3];

    Smooth3D::Real3 axe1 = (1. - lambda) * (z1 - z0);
    Smooth3D::Real3 dir = math::vecMul(axe1, normal);
    new_nodes[(i + 1) % 3] = z0 + axe1 + dir * t;
  }
  // Deuxieme passe
  for (int i = 0; i < 3; i++) {
    Smooth3D::Real3 z0 = new_nodes[i];
    Smooth3D::Real3 z1 = new_nodes[(i + 1) % 3];

    Smooth3D::Real3 axe1 = (1. - lambda) * (z0 - z1);
    Smooth3D::Real3 dir = -math::vecMul(axe1, normal);

    compute_nodes[i] = z1 + axe1 + dir * t;
  }

  // Scaling !
  Smooth3D::Real perimeter_update = (compute_nodes[1] - compute_nodes[0]).abs()
      + (compute_nodes[2] - compute_nodes[1]).abs()
      + (compute_nodes[2] - compute_nodes[0]).abs();
  Smooth3D::Real scale = perimeter / perimeter_update;

  compute_nodes[0] = centroid + (compute_nodes[0] - centroid) * scale;
  compute_nodes[1] = centroid + (compute_nodes[1] - centroid) * scale;
  compute_nodes[2] = centroid + (compute_nodes[2] - centroid) * scale;

  return compute_nodes;
}

static FourNodes computeGetMeQuadrangle3D(gmds::Face& face,
    Smooth3D::Real theta, Smooth3D::Real lambda) {
  FourNodes compute_nodes, new_nodes;

  std::vector<gmds::Node> quadrangles_nodes = face.get<gmds::Node>();
  for (int i = 0; i < 4; i++) {
    gmds::Node current_node = quadrangles_nodes[i];
    compute_nodes[i] = Smooth3D::Real3(current_node.X(), current_node.Y(),
	current_node.Z());
    // We map the nodes in 2D, Z0 => (0.0)

  }
  Smooth3D::Real3 centroid = (compute_nodes[0] + compute_nodes[1]
      + compute_nodes[2] + compute_nodes[3]) * (1.0 / 4.0);
  Smooth3D::Real perimeter = (compute_nodes[1] - compute_nodes[0]).abs()
      + (compute_nodes[2] - compute_nodes[1]).abs()
      + (compute_nodes[3] - compute_nodes[2]).abs()
      + (compute_nodes[0] - compute_nodes[3]).abs();
  Smooth3D::Real3 normal = math::vecMul(compute_nodes[2] - compute_nodes[0],
      compute_nodes[3] - compute_nodes[1]);
  normal = normal / normal.abs();

  Smooth3D::Real t = std::tan(theta);
  // premiere passe

  for (int i = 0; i < 4; i++) {
    Smooth3D::Real3 z0 = compute_nodes[i];
    Smooth3D::Real3 z1 = compute_nodes[(i + 1) % 4];

    Smooth3D::Real3 axe1 = (1. - lambda) * (z1 - z0);
    Smooth3D::Real3 dir = math::vecMul(axe1, normal);

    new_nodes[(i + 1) % 4] = z0 + axe1 + dir * t;

  }

  // Deuxieme passe
  for (int i = 0; i < 4; i++) {
    Smooth3D::Real3 z0 = new_nodes[i];
    Smooth3D::Real3 z1 = new_nodes[(i + 1) % 4];

    Smooth3D::Real3 axe1 = (1. - lambda) * (z0 - z1);
    Smooth3D::Real3 dir = -math::vecMul(axe1, normal);

    compute_nodes[i] = z1 + axe1 + dir * t;
  }

  // Scaling !
  Smooth3D::Real perimeter_update = (compute_nodes[1] - compute_nodes[0]).abs()
      + (compute_nodes[2] - compute_nodes[1]).abs()
      + (compute_nodes[3] - compute_nodes[2]).abs()
      + (compute_nodes[0] - compute_nodes[3]).abs();
  Smooth3D::Real scale = perimeter / perimeter_update;

  compute_nodes[0] = centroid + (compute_nodes[0] - centroid) * scale;
  compute_nodes[1] = centroid + (compute_nodes[1] - centroid) * scale;
  compute_nodes[2] = centroid + (compute_nodes[2] - centroid) * scale;
  compute_nodes[3] = centroid + (compute_nodes[3] - centroid) * scale;

  return compute_nodes;
}

static void computeN2N(gmds::Mesh& mesh) {

  std::map<gmds::TCellID, std::set<gmds::TCellID> > n2n;

  gmds::Mesh::face_iterator it_faces = mesh.faces_begin();

  while (!it_faces.isDone()) {
    gmds::Face f = it_faces.value();

    std::vector<gmds::Node> f_nodes = f.get<gmds::Node>();
    for (unsigned int i = 0; i < f_nodes.size(); i++) {
      gmds::Node ni = f_nodes[i];
      gmds::Node nj, nk;
      f.getAdjacentNodes(ni, nj, nk);
      n2n[ni.getID()].insert(nj.getID());
      n2n[ni.getID()].insert(nk.getID());
    }
    it_faces.next();
  }

  std::map<gmds::TCellID, std::set<gmds::TCellID> >::iterator it_map;
  for (it_map = n2n.begin(); it_map != n2n.end(); it_map++) {
    gmds::Node current_node = mesh.get<gmds::Node>(it_map->first);
    std::set<gmds::TCellID> set_adj_nodes = it_map->second;
    std::vector<gmds::TCellID> adj_nodes;
    adj_nodes.insert(adj_nodes.end(), set_adj_nodes.begin(),
	set_adj_nodes.end());
    current_node.set<gmds::Node>(adj_nodes);
  }
}
*/

extern "C" int S3_GETMe2D(const double alpha, const double beta,
    int_type nb_cells, int_type nb_nodes, const int_type * nb_nodes_per_cell,
    const int_type *nodes_number, double *x, double *y, double *z,
    const double * weights, const double * relax, int_type n_iter) {
    std::cerr<<"S3_GETMe2D is no longer available."<<std::endl;
    return -1;
  /*
  const int kTGMDSMask = gmds::DIM3 | gmds::N | gmds::F | gmds::F2N | gmds::N2N;
  gmds::Mesh internal_mesh(kTGMDSMask);

  gmds::Node* temp_array_nodes = new gmds::Node[nb_nodes];

  for (int vtx = 0; vtx < nb_nodes; ++vtx) {
    temp_array_nodes[vtx] = gmds::Node();
  }
  const int_type * ptr_connectivity = nodes_number;
  for (int cell = 0; cell < nb_cells; ++cell) {
    std::vector<gmds::Node> nodes_cell;
    for (int vtx = 0; vtx < nb_nodes_per_cell[cell]; ++vtx) {
      int which_node = *ptr_connectivity++;
      if (temp_array_nodes[which_node].getID() == gmds::NullID) {
	temp_array_nodes[which_node] = internal_mesh.newNode(x[which_node],
	    y[which_node], z[which_node]);

      }
      nodes_cell.push_back(temp_array_nodes[which_node]);
    }
    internal_mesh.newFace(nodes_cell);
  }

  gmds::MeshDoctor mesh_doc(&internal_mesh);
  mesh_doc.updateUpwardConnectivity();
  computeN2N(internal_mesh);

  Smooth3D::Real theta_max = (1.0 / 6.0) * M_PI;
  Smooth3D::Real3 *new_pos = new Smooth3D::Real3[internal_mesh.getNbNodes()];
  Smooth3D::Real *weights_pos = new Smooth3D::Real[internal_mesh.getNbNodes()];
  Smooth3D::Real3 *bi = new Smooth3D::Real3[internal_mesh.getNbNodes()];
  Smooth3D::Real3 *orig = new Smooth3D::Real3[internal_mesh.getNbNodes()];

  for (int local_id = 0; local_id < internal_mesh.getNbNodes(); local_id++) {
    bi[local_id] = Smooth3D::Real3::null();
  }
  for (int vtx = 0; vtx < nb_nodes; vtx++) {
    if (temp_array_nodes[vtx].getID() != gmds::NullID && relax[vtx] != 0.0) {
      gmds::Node current_node = temp_array_nodes[vtx];
      int local_id = current_node.getID();
      orig[local_id] = Smooth3D::Real3(current_node.X(), current_node.Y(),
	  current_node.Z());
    }
  }
  for (int iter_num = 0; iter_num < n_iter; iter_num++) {

    for (int vtx = 0; vtx < internal_mesh.getNbNodes(); vtx++) {
      new_pos[vtx] = Smooth3D::Real3::null();
      weights_pos[vtx] = 0.0;
    }

    for (gmds::Mesh::face_iterator it = internal_mesh.faces_begin();
	!it.isDone(); it.next()) {
      gmds::Face current_face = it.value();

      FourNodes new_position;

      Smooth3D::Real quality = computeMeanRatio(current_face);
      Smooth3D::Real weight = 0.1 + (1.0 - quality);
      Smooth3D::Real theta = (1.0 - quality) * theta_max;
      Smooth3D::Real lambda, lambda_max;
      if (current_face.getNbNodes() == 3)
	lambda_max = 0;
      else
	lambda_max = 0.25;
      lambda = lambda_max * (1.0 - quality);

      if (current_face.getNbNodes() == 3) {
	std::vector<gmds::Node> triangle_nodes = current_face.get<gmds::Node>();
	new_position = computeGetMeTriangle3D(current_face, theta, lambda);
	for (int i = 0; i < 3; i++) {
	  new_pos[triangle_nodes[i].getID()] += weight * new_position[i];
	  weights_pos[triangle_nodes[i].getID()] += weight;
	}
      } else if (current_face.getNbNodes() == 4) {
	std::vector<gmds::Node> quadrangle_nodes =
	    current_face.get<gmds::Node>();
	new_position = computeGetMeQuadrangle3D(current_face, theta, lambda);
	for (int i = 0; i < 4; i++) {
	  new_pos[quadrangle_nodes[i].getID()] += weight * new_position[i];
	  weights_pos[quadrangle_nodes[i].getID()] += weight;
	}
      }
    }

    // We use the HC algorithm from "Improved Laplacian Smoothing of Noisy Surface meshes"
    // Eurographics '99, J. Vollmer, R. Mencl and H. Muller

    // HC First pass, compute the new pos and the b !
    for (int vtx = 0; vtx < nb_nodes; vtx++) {
      gmds::Node current_node = temp_array_nodes[vtx];
      int local_id = current_node.getID();
      if (local_id != gmds::NullID && relax[vtx] != 0.0) {
	Smooth3D::Real3 q = Smooth3D::Real3(current_node.X(), current_node.Y(),
	    current_node.Z());
	Smooth3D::Real3 pos = q;
	if (weights_pos[local_id] != 0.0) {
	  pos = new_pos[local_id] / weights_pos[local_id];
	  current_node.setXYZ(pos.x(), pos.y(), pos.z());
	}
      bi[local_id] = pos - (alpha * orig[local_id] + (1.0- alpha) * q);
    }
  }

  // HC second pass, ajust the pos regarding the b of the neighbourhood
  for (int vtx = 0; vtx < nb_nodes; vtx++) {
    gmds::Node current_node = temp_array_nodes[vtx];
    int local_id = current_node.getID();
    if (local_id != gmds::NullID && relax[vtx] != 0.0) {

      Smooth3D::Real3 somme = Smooth3D::Real3::null();
      Smooth3D::Real number = 0.;

      std::vector<gmds::Node> neighbor_nodes = current_node.get<gmds::Node>();
      for (int vtx2 = 0; vtx2 < neighbor_nodes.size(); ++vtx2) {
	gmds::Node voisin = neighbor_nodes[vtx2];
	int voisin_id = voisin.getID();
	somme = somme + bi[voisin_id];
	number += 1.0;
      }

      if (number != 0.) {
	Smooth3D::Real3 pos = Smooth3D::Real3(current_node.X(),
	    current_node.Y(), current_node.Z());
	pos = pos - (beta * bi[local_id] + (1.0 - beta) / number * somme);
	current_node.setXYZ(pos.x(), pos.y(), pos.z());
      }
    }
  }
}

for (int vtx = 0; vtx < nb_nodes; vtx++) {
  if (temp_array_nodes[vtx].getID() != gmds::NullID) {
    gmds::Node current_node = temp_array_nodes[vtx];
    double rel = relax[vtx];

    if (rel != 0.0) {
      double xo = current_node.X();
      double yo = current_node.Y();
      double zo = current_node.Z();

      x[vtx] = (1.0 - rel) * x[vtx] + rel * xo;
      y[vtx] = (1.0 - rel) * y[vtx] + rel * yo;
      z[vtx] = (1.0 - rel) * z[vtx] + rel * zo;
    }
  }
}

delete[] orig;
delete[] bi;
delete[] weights_pos;
delete[] new_pos;

delete[] temp_array_nodes;

return (0);
*/
}

/*
void testGetMeUnitTri(void) {
Smooth3D::Real theta = M_PI / 6.0;
Smooth3D::Real lambda = 0.0;

const int kTGmdsMask = gmds::DIM3 | gmds::N | gmds::F | gmds::F2N;
gmds::Mesh internal_mesh(kTGmdsMask);

// we start with a square triangle
gmds::Node noeuds[3];
noeuds[0] = internal_mesh.newNode(0.0, 0.0, 0.0);
noeuds[1] = internal_mesh.newNode(3.0, 0.0, 0.0);
noeuds[2] = internal_mesh.newNode(0.0, 4.0, 0.0);

Smooth3D::Real xg = 3.0 / 3.0, yg = 4.0 / 3.0, zg = 0.0;

gmds::Face face = internal_mesh.newTriangle(noeuds[0], noeuds[1], noeuds[2]);
for (int i = 0; i < 10; i++) {

  Smooth3D::Real qual = computeMeanRatio(face);
  std::cout << "iter " << i << " qual= " << qual << " Centroid=(" << xg << ","
      << yg << "," << zg << ")" << std::endl;
  FourNodes newpos = computeGetMeTriangle3D(face, theta, lambda);

  xg = yg = zg = 0.0;

  for (int a = 0; a < 3; a++) {
    xg += newpos[a].x();
    yg += newpos[a].y();
    zg += newpos[a].z();
    noeuds[a].setXYZ(newpos[a].x(), newpos[a].y(), newpos[a].z());
  }
  xg = xg / 3.0;
  yg = yg / 3.0;
  zg = zg / 3.0;
}

}

void testGetMeUnitQuad(void) {
Smooth3D::Real theta = M_PI / 6.0;
Smooth3D::Real lambda = 0.25;

const int kTGMDSMask = gmds::DIM3 | gmds::N | gmds::F | gmds::F2N;
gmds::Mesh internal_mesh(kTGMDSMask);

// we start with a square triangle
gmds::Node noeuds[4];
noeuds[0] = internal_mesh.newNode(0.0, 0.0, 0.0);
noeuds[1] = internal_mesh.newNode(3.0, 0.0, 0.0);
noeuds[2] = internal_mesh.newNode(0.0, 4.0, 0.0);
noeuds[3] = internal_mesh.newNode(-1.0, 1.0, 0.0);

Smooth3D::Real xg = 2.0 / 4.0, yg = 5.0 / 4.0, zg = 0.0;

gmds::Face face = internal_mesh.newQuad(noeuds[0], noeuds[1], noeuds[2],
    noeuds[3]);
for (int i = 0; i < 10; i++) {

  Smooth3D::Real qual = computeMeanRatio(face);
  std::cout << "iter " << i << " qual= " << qual << " Centroid=(" << xg << ","
      << yg << "," << zg << ")" << std::endl;
  FourNodes newpos = computeGetMeQuadrangle3D(face, theta, lambda);

  xg = yg = zg = 0.0;

  for (int a = 0; a < 4; a++) {
    xg += newpos[a].x();
    yg += newpos[a].y();
    zg += newpos[a].z();
    noeuds[a].setXYZ(newpos[a].x(), newpos[a].y(), newpos[a].z());
  }
  xg = xg / 4.0;
  yg = yg / 4.0;
  zg = zg / 4.0;
}

}

void testGetMeUnitQuad2(void) {
Smooth3D::Real theta = M_PI / 6.0;
Smooth3D::Real lambda = 0.25;

const int kTGMDSMask = gmds::DIM3 | gmds::N | gmds::F | gmds::F2N;
gmds::Mesh internal_mesh(kTGMDSMask);

// we start with a non convex quadrangle

gmds::Node noeuds[4];
noeuds[0] = internal_mesh.newNode(0.0, 0.0, 0.0);
noeuds[1] = internal_mesh.newNode(3.0, 0.0, 0.0);
noeuds[2] = internal_mesh.newNode(0.0, 4.0, 0.0);
noeuds[3] = internal_mesh.newNode(1.0, 1.0, 0.0);

Smooth3D::Real xg = 4.0 / 4.0, yg = 5.0 / 4.0, zg = 0.0;

gmds::Face face = internal_mesh.newQuad(noeuds[0], noeuds[1], noeuds[2],
    noeuds[3]);
for (int i = 0; i < 10; i++) {

  Smooth3D::Real qual = computeMeanRatio(face);
  std::cout << "iter " << i << " qual= " << qual << " Centroid=(" << xg << ","
      << yg << "," << zg << ")" << std::endl;

  FourNodes newpos = computeGetMeQuadrangle3D(face, theta, lambda);

  xg = yg = zg = 0.0;

  for (int a = 0; a < 4; a++) {
    xg += newpos[a].x();
    yg += newpos[a].y();
    zg += newpos[a].z();
    noeuds[a].setXYZ(newpos[a].x(), newpos[a].y(), newpos[a].z());
  }
  xg = xg / 4.0;
  yg = yg / 4.0;
  zg = zg / 4.0;
}

}
*/