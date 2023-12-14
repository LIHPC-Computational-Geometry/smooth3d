/*
 * s3_conditionnumber2D.cpp
 *
 *  Created on: 13 nov. 2013
 *      Author: weilljc
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cstdlib>		// pour rand

#include "smooth3D/smooth.h"
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>

#include "math/Real3.h"

#include "math/BFGS.h"
#include "math/ConjugateGradientPr.h"
#include "ObjectiveCondition2D.h"
#include "ObjectiveLag2D.h"

extern "C" int S3_ConditionNumber2D(
		int_type nb_cells, int_type nb_nodes,
		const int_type * nb_nodes_per_cell,
		const int_type *nodes_number,
		double *x, double *y, double *z,
		const double * weights, const double * relax,
		int_type n_iter)
{
  const int tGMDSMask = gmds::DIM3 | gmds::N | gmds::F | gmds::F2N | gmds::N2F;
  gmds::Mesh  internal_mesh(tGMDSMask);

  gmds::Node * temp_array_nodes = new gmds::Node [nb_nodes];

  for (int vtx = 0; vtx < nb_nodes; ++vtx) {
    temp_array_nodes[vtx] = internal_mesh.newNode(x[vtx], y[vtx], z[vtx]);
  }
  const int_type * ptr_connectivity = nodes_number;
  for (int cell = 0; cell < nb_cells; ++cell) {
    std::vector<gmds::Node> nodes_cell;
    for (int vtx = 0; vtx < nb_nodes_per_cell[cell]; ++vtx)
      nodes_cell.push_back(temp_array_nodes[*ptr_connectivity++]);
    internal_mesh.newFace(nodes_cell);
  }

  gmds::MeshDoctor mesh_doc(&internal_mesh);
  mesh_doc.updateUpwardConnectivity();

  delete[] temp_array_nodes;

  Smooth3D::ObjectiveCondition2D condition;

  Smooth3D::Real3 pos[4]; // 4 nodes per faces at most ! Beware

  Smooth3D::Real max_f = 0.0;
  // LOOP 1 : FOR ALL ITER
  for (int iter = 0; iter < n_iter; iter++) {

    // LOOP 2 : FOR ALL NODES
    for (auto i : internal_mesh.nodes()){
      auto current_node = internal_mesh.get<gmds::Node>(i);
      if (relax[current_node.id()] != 0.0) {

	// We are only interested on Relaxed nodes
	condition.set_node(current_node);
	int nb_faces = current_node.nbFaces();

	// First face connected to the node is chosen ramdomly
	int face_int = std::rand() % nb_faces;
	std::vector<gmds::Face> allFaces = current_node.get<gmds::Face>();
	bool valide = false;

	// LOOP 3 : FOR ALL FACES CONNECTED TO THE NODE
	for (int face_i = 0; !valide && face_i < nb_faces; face_i++) {
	  gmds::Face current_face = allFaces[(face_i + face_int) % nb_faces];
	  std::vector<gmds::Node > all_nodes = current_face.get<gmds::Node>();
	  int nb_vtx = all_nodes.size();
	  if (nb_vtx == 3 || nb_vtx == 4) {
	    // We are having here either a triangle or a quadrangle
	    for (int i = 0; i < nb_vtx; ++i)
	      if (all_nodes[i] == current_node) {
		for (int j = 0; j < nb_vtx; ++j) {
		  gmds::Node tempNode = all_nodes[(i + j) % nb_vtx];
		  pos[j] = Smooth3D::Real3(tempNode.X(), tempNode.Y(), tempNode.Z());
		}
		break;
	      }

	    Smooth3D::LagrangeMapping2D * p = 0;
	    if (nb_vtx == 3)
	    	p = new Smooth3D::TriangleLagrangeMapping2D(pos[0], pos[1],
	    			pos[2]);
	    else
	    	p = new Smooth3D::QuadrangleLagrangeMapping2D(pos[0], pos[1],
	    			pos[2], pos[3]);

	    Smooth3D::ObjectiveLag2D my_objective(&condition, p);
	    Optimize::BFGS<Smooth3D::Real2> algo_opt(&my_objective, 2);
	    Smooth3D::Real2 x(0.0, 0.0);
	    Smooth3D::Real step_size = 0.01;

	    algo_opt.globalSet(x, step_size, 0.01);

	    Smooth3D::Real depart = algo_opt.value();
	    Smooth3D::Real best_val = depart;
	    Smooth3D::Real2 best = x;

	    for (int iter = 0; iter < 10; iter++) {
	    	int status = algo_opt.globalIterate();
	    	x = algo_opt.solution();
	    	if (p->isvalid(x) && algo_opt.value() < best_val) {
	    		best_val = algo_opt.value();
	    		best = x;
	    	}
	    	if (status != Optimize::OPT_SUCCESS
	    			|| algo_opt.convergenceP(1.0e-4)) {
	    		break;
	    	}

	    }

	    x = best;

	    if (p->isvalid(x)) {
	    	Smooth3D::Real3 real_pos = p->eval(x);
	    	Smooth3D::Real val = best_val;

	    	if (val > max_f)
	    		max_f = val;

	    	if (val < depart) {
	    		current_node.setXYZ(real_pos.x(), real_pos.y(), real_pos.z());
	    		valide = true;
	    	}
	    }

	    delete p;
	  }
	}
      }
    }
  }


  for (auto i : internal_mesh.nodes()){
       auto current_node = internal_mesh.get<gmds::Node>(i);
       int vtx = current_node.id();
       double rel = relax[vtx];

       if (rel  != 0.0)
       {
	 double xo = current_node.X();
	 double yo = current_node.Y();
	 double zo = current_node.Z();


	  x[vtx] = (1.0 - rel) * x[vtx] + rel * xo;
	  y[vtx] = (1.0 - rel) * y[vtx] + rel * yo;
	  z[vtx] = (1.0 - rel) * z[vtx] + rel * zo;
       }

  }
  return 1;
}

