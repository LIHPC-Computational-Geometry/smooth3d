/*
 * This will be the header.,...

 */

#ifndef SMOOTH3D
#define SMOOTH3D

#include "machine_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int S3_laplace(int_type nb_cells, int_type nb_nodes,
    const int_type *nb_nodes_per_cell, const int_type *nodes_number, double *x,
    double *y, double *z, const double *weights, const double *relax,
    int_type n_iter);
int S3_Tipton(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double * x, double *y, double * z, const double * weights,
    const double * relax, int_type n_iter);
int S3_Jun(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double * x, double *y, double * z, const double * weights,
    const double * relax, int_type n_iter);
int S3_ConditionNumber(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter);

int S3_InverseMeanRatio(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter);

int S3_ConditionNumber2D(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type *nodes_number, double *x,
    double *y, double *z, const double * weights, const double * relax,
    int_type n_iter);

// S3_GETMe2D is no longer available.
int S3_GETMe2D(const double alpha, const double beta, int_type nb_cells,
    int_type nb_nodes, const int_type * nb_nodes_per_cell,
    const int_type *nodes_number, double *x, double *y, double *z,
    const double * weights, const double * relax, int_type n_iter);

// S3_GETMe3D is no longer available. 
// - fixed_nodes : array of size nb_nodes carrying 0 or 1
// - n_iterSimult : number of iterations for the simultaneous GETMe
// - n_iterSeq : number of iterations for the sequential GETMe
int S3_GETMe3D(const int_type nb_cells,
    int_type nb_nodes, const int_type * nb_nodes_per_cell,
    const int_type *nodes_indices, double *x, double *y, double *z,
    const int_type *fixed_nodes, double quality_threshold, int_type n_iterSimult, int_type n_iterSeq);

int S3_Orthogonal2D(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter);

int S3_Jun2D(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double * x, double *y, double * z, const double * weights,
    const double * relax, int_type n_iter);

#ifdef __cplusplus
}
;
#endif

#endif
