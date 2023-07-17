/*
 * Generic 3D Geometry to be used with the NPROJ project
 *
 * All you wanted to know about Volume / Centroids
 *                                / Decomposition in TetraHedras
 * $Id$
 */

#include "smooth3D/smooth.h"



/* Return the Faces Array of a given Cell */

const int kS3NbFacesCell[9] =
{
  0, 0, 0, 0, 4, 5, 5, 0, 6};

int 
s3Faces3DCell (s3_cell_3d_t const *cell_3d, s3_face_t faces[6])
{
#define face_copy_vect(n, i, j) {      \
      faces[n].ind_node [i]= j;\
}


#define face_triangle(n, a, b, c) {   \
       face_copy_vect (n, 0, a); \
       face_copy_vect (n, 1, b); \
       face_copy_vect (n, 2, c); \
       faces[n].nb_nodes = FACE_TRI;    \
}

#define face_quadrangle(n, a, b, c, d) {\
        face_copy_vect (n, 0, a); \
	face_copy_vect (n, 1, b); \
        face_copy_vect (n, 2, c); \
        face_copy_vect (n, 3, d); \
        faces[n].nb_nodes = FACE_QUA;    \
}

  switch (cell_3d->nb_nodes)
    {
    case R3D_TETRA:
      face_triangle (0, 0, 1, 3);
      face_triangle (1, 0, 3, 2);
      face_triangle (2, 0, 2, 1);
      face_triangle (3, 1, 2, 3);
      return 4;
    case R3D_PYRAMID:
      face_triangle (0, 0, 1, 4);
      face_triangle (1, 0, 4, 3);
      face_triangle (2, 1, 2, 4);
      face_triangle (3, 3, 4, 2);
      face_quadrangle (4, 0, 3, 2, 1);
      return 5;
    case R3D_WEDGE:
      /*                        
       *                       
       *                        
       *                      3 
       *                       /|\
       *                      / | \
       *                     /  |  \
       *                  4 |-------| 5
       *                    |   |   |
       *                    |  0|   |
       *                    |  / \  |
       *                    | /   \ |
       *                    |/     \|
       *                  1 !-------| 2
       *                                   
       *                   
       *                 
       */
      face_triangle (0, 0, 2, 1);
      face_triangle (1, 3, 4, 5);
      face_quadrangle (2, 0, 3, 5, 2);
      face_quadrangle (3, 0, 1, 4, 3);
      face_quadrangle (4, 1, 2, 5, 4);
      return 5;
    case R3D_HEXA:
      /*                         
       *                        
       *                         
       *                       4 |--------------|7
       *                        /|             /|
       *                       / |            / |
       *                      /  |         6 /  |
       *                   5 |--------------|   |
       *                     |   |          |   |
       *                     |   |----------|---|
       *                     |  / 0         |  / 3        
       *                     | /            | /
       *                     |/             |/
       *                   1 !--------------| 2
       *                                    
       *                  
       */

      face_quadrangle (0, 0, 1, 5, 4);
      face_quadrangle (1, 0, 3, 2, 1);
      face_quadrangle (2, 0, 4, 7, 3);
      face_quadrangle (3, 1, 2, 6, 5);
      face_quadrangle (4, 3, 7, 6, 2);
      face_quadrangle (5, 4, 5, 6, 7);
      return 6;
    default:
      fprintf (stderr, "Internal Error in the Smoothing routines\n");
      assert (1 == 0);
    }
  return -1;
}


int 
s3Edges3DCell (s3_cell_3d_t const *cell_3d, s3_edge_t edges[12])
{
#define edge_copy_vect(n, f, t) {      \
      edges[n].from = f, edges[n].to = t;\
}


  switch (cell_3d->nb_nodes)
    {
    case R3D_TETRA:
      edge_copy_vect (0, 0, 1);
      edge_copy_vect (1, 1, 2);
      edge_copy_vect (2, 2, 0);
      edge_copy_vect (3, 0, 3);
      edge_copy_vect (4, 1, 3);
      edge_copy_vect (5, 3, 2);
      return 6;
    case R3D_PYRAMID:
      edge_copy_vect (0, 0, 1);
      edge_copy_vect (1, 1, 2);
      edge_copy_vect (2, 2, 3);
      edge_copy_vect (3, 3, 0);
      edge_copy_vect (4, 0, 4);
      edge_copy_vect (5, 1, 4);
      edge_copy_vect (6, 2, 4);
      edge_copy_vect (7, 3, 4);
      return 8;
    case R3D_WEDGE:

      /*                        
       *                       
       *                        
       *                      3 
       *                       /|\
       *                      / | \
       *                     /  |  \
       *                  4 |-------| 5
       *                    |   |   |
       *                    |  0|   |
       *                    |  / \  |
       *                    | /   \ |
       *                    |/     \|
       *                  1 !-------| 2
       *                                   
       *                   
       *                 
       */

      edge_copy_vect (0, 0, 1);
      edge_copy_vect (1, 1, 2);
      edge_copy_vect (2, 2, 0);
      edge_copy_vect (3, 3, 4);
      edge_copy_vect (4, 4, 5);
      edge_copy_vect (5, 5, 3);
      edge_copy_vect (6, 1, 4);
      edge_copy_vect (7, 2, 5);
      edge_copy_vect (8, 0, 3);
      return 9;

    case R3D_HEXA:
      /*                         
       *                        
       *                         
       *                       4 |--------------|7
       *                        /|             /|
       *                       / |            / |
       *                      /  |         6 /  |
       *                   5 |--------------|   |
       *                     |   |          |   |
       *                     |   |----------|---|
       *                     |  / 0         |  / 3        
       *                     | /            | /
       *                     |/             |/
       *                   1 !--------------| 2
       *                                    
       *                  
       */  
      edge_copy_vect (0, 0, 1);
      edge_copy_vect (1, 1, 2);
      edge_copy_vect (2, 2, 3);
      edge_copy_vect (3, 3, 0);
      edge_copy_vect (4, 4, 5);
      edge_copy_vect (5, 5, 6);
      edge_copy_vect (6, 6, 7);
      edge_copy_vect (7, 7, 4);
      edge_copy_vect (8, 0, 4);
      edge_copy_vect (9, 1, 5);
      edge_copy_vect (10, 2, 6);
      edge_copy_vect (11, 3, 7);
      return 12;
    default:
      fprintf (stderr, "Internal Error in the Smoothing routines\n");
      assert (1 == 0);
    }
  return -1;
}
