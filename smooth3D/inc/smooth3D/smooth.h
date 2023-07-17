/*
 *
 * $Id$
 * Version originale des fonctions helpers pour l'interface C
 */

#ifndef SMOOTH_INCLUDE_H
#define SMOOTH_INCLUDE_H

#include <stdio.h>
#include <fcntl.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <values.h>

#include <machine_types.h>
#include "smooth3D/smooth3d.h"

/*
 * Definition of 3d Objects
 */

enum S3_GEOM {
  R3D_TETRA = 4,
  R3D_PYRAMID = 5,
  R3D_WEDGE = 6,
  R3D_HEXA = 8,
  R3D_UNKNOWN = -1
};

typedef enum S3_GEOM s3_geom_3d_t;

enum S3FACES {
  FACE_TRI = 3,
  FACE_QUA = 4
};
typedef enum S3FACES S3face_kind_t;

/* Now a 3d Cell */
struct S3_CELL_3D {
  s3_geom_3d_t nb_nodes;
  int_type vtx[8];
};

typedef struct S3_CELL_3D s3_cell_3d_t;

/* a 3d Face... Can be a triangle (planar of course) or a quad which is
 supposely not always planar
 */

struct S3_FACE {
  S3face_kind_t nb_nodes;
  int_type ind_node[4]; /* Indices in the Cell */
};

typedef struct S3_FACE s3_face_t;

struct S3_EDGE {
  int_type from, to;
};

typedef struct S3_EDGE s3_edge_t;

/* src/s3_3dgeom.c */
int s3Faces3DCell(s3_cell_3d_t const *cell_3d, s3_face_t faces[6]);
int s3Edges3DCell(s3_cell_3d_t const *cell_3d, s3_edge_t edges[12]);
/* src/s3mem.c */
void s3Free(void *p);
void *s3Malloc(size_t taille);
void *s3Realloc(void *p1, size_t taille);

#endif
