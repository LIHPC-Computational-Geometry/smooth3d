/*
 * For the projection Project
 *
 * Memory access... So can be changed easily if needed
 */

#include "smooth3D/smooth.h"

void s3Free (void * p)
{
#ifdef DEBUG_MEM
  fprintf(stderr, "free %x\n", (int) p);
#endif
  free(p);
}

void * s3Malloc (size_t taille)
{
 void * p;
 p = malloc (taille);
#ifdef DEBUG_MEM
   fprintf(stderr, "Alloc size = %d => %x\n", (int) taille, (int) p);
#endif
  return p;
}

void * s3Realloc (void * p1, size_t taille)
{
  void * p = realloc(p1, taille);
#ifdef DEBUG_MEM
   fprintf(stderr, "Realloc %x size = %d => %x\n",
	   (int) p1, (int) taille, (int) p);
#endif
  return p;
}
