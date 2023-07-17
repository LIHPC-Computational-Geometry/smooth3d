/*
 * Smooth3DMesh.h
 *
 *  Created on: 21 févr. 2013
 *      Author: weilljc
 */

#ifndef SMOOTH3DMESH_H_
#define SMOOTH3DMESH_H_

#include <MeshImpl.hpp>
#include <machine_types.h>


class Smooth3DMesh : public Mesquite2::MeshImpl{
public:
  Smooth3DMesh();
  virtual ~Smooth3DMesh();

  void initSmooth3DMesh(int_type nb_cells, int_type nb_nodes,
	   const int_type * nb_nodes_per_cell,
	   const int_type * nodes_number,
	   double *x, double *y, double *z);

  void epilogue(double *x, double *y, double *z, const double * relax);

  /* Les noeuds à relax == 0.0 sont mis à fixed */

  void markNonRelaxedFixed(const double *relax);

};


#endif /* SMOOTH3DMESH_H_ */
