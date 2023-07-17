/*
 * Methode du Laplacien
 * $Id$
 * $Log$
 * Revision 1.2  2000/11/16 11:47:49  weilljc
 * Version 0.0.3
 *
 *
 *    Ajout de S3_Tipton...
 *
 * Revision 1.1.1.1  2000/10/25 12:28:39  weilljc
 * Lissage de maillages tri-dim
 *
 */


#include "smooth3D/smooth.h"

int S3_laplace (int_type nb_cells, int_type nb_nodes,
		const int_type  * nb_nodes_per_cell,
		const int_type  * nodes_number,
		double * x, double *y, double * z,
		const double  * weights, const double  * relax,
		int_type n_iter)
{
  
  
  
  int err = 0;			/*  Error Flags... 0 means error, 1 means oK*/
  double * sum_x = NULL, * sum_y =NULL, * sum_z=NULL;
  double * total_weight;
  int_type vtx, cell;
  const int_type  * f_node_cell;
  int iter;

  total_weight = (double *) s3Malloc (nb_nodes * sizeof(total_weight[0]));
  if (total_weight == NULL)
    goto epilogue;
  sum_x = (double *) s3Malloc (nb_nodes * sizeof(sum_x[0]));
  if (sum_x == NULL)
    goto epilogue;
  sum_y = (double *) s3Malloc (nb_nodes * sizeof(sum_y[0]));
  if (sum_y == NULL)
    goto epilogue;
  sum_z = (double *) s3Malloc (nb_nodes * sizeof(sum_z[0]));
  if (sum_z == NULL)
    goto epilogue;
  
  for (vtx = 0; vtx < nb_nodes; ++vtx)
    {
      total_weight[vtx] = 0.0;
      sum_x[vtx] = sum_y[vtx] = sum_z[vtx] = 0.0;
    }

  for (iter = 0; iter < n_iter; ++iter) 
    {
      for (f_node_cell = nodes_number, cell = 0; cell < nb_cells; ++cell)
	{
	  s3_cell_3d_t my_cell;
	  s3_edge_t my_edges[12];
	  int nb_edges, edge;

	  my_cell.nb_nodes = nb_nodes_per_cell[cell];
	  for (vtx = 0; vtx < nb_nodes_per_cell[cell]; ++vtx)
	    my_cell.vtx[vtx] = *(f_node_cell++);

	  nb_edges = s3Edges3DCell(&my_cell, my_edges);
	  if (nb_edges == -1)
	    goto epilogue;


	  for (edge = 0; edge < nb_edges; ++edge)
	    {
	      int vtx = my_edges[edge].from;
	      int vtx1 = my_edges[edge].to;
	      int id1 = my_cell.vtx[vtx];
	      int id2 =  my_cell.vtx[vtx1];
	      sum_x[id2] += weights[id1] * x[id1];
	      sum_y[id2] += weights[id1] * y[id1];
	      sum_z[id2] += weights[id1] * z[id1];
		  
	      sum_x[id1] += weights[id2] * x[id2];
	      sum_y[id1] += weights[id2] * y[id2];
	      sum_z[id1] += weights[id2] * z[id2];
	      total_weight[id2] += weights[id1];
	      total_weight[id1] += weights[id2];
	    }
	}
      for (vtx = 0; vtx < nb_nodes; ++vtx)
	{
	  double rel = relax[vtx];
	  if (rel != 0.0 && total_weight[vtx] != 0)
	    {
	      x[vtx] = (1-rel)* x[vtx] + rel * (sum_x[vtx] / total_weight[vtx]);
	      y[vtx] = (1-rel)* y[vtx] + rel * (sum_y[vtx] / total_weight[vtx]);
	      z[vtx] = (1-rel)* z[vtx] + rel  * (sum_z[vtx] / total_weight[vtx]);
	    }
	}
    }

	

  err = 1;			/* All went well */

  epilogue:
  if (sum_x != NULL)
	  s3Free(sum_x);
  if (sum_y != NULL)
	  s3Free(sum_y);
  if (sum_z != NULL)
	  s3Free(sum_z);
  if (total_weight != NULL)
	  s3Free(total_weight);

  return err;
}
