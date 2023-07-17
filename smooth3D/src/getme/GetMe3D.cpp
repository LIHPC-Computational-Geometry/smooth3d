/*
 * GetMe3D.cpp
 */

#include "smooth3D/smooth.h"

#include <iostream>

extern "C" int S3_GETMe3D(const int_type nb_cells,
	    int_type nb_nodes, const int_type * nb_nodes_per_cell,
	    const int_type *nodes_indices, double *x, double *y, double *z,
	    const int_type *fixed_nodes, double quality_threshold, int_type n_iterSimult, int_type n_iterSeq)
{
	std::cerr<<"S3_GETMe3D is no longer available."<<std::endl;
        return -1;
}
