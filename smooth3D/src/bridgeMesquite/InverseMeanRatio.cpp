/*
 * InverseMeanRation.cpp
 *
 *  Created on: 22 févr. 2013
 *      Author: weilljc
 */

#include "Smooth3DMesh.h"
#include <InstructionQueue.hpp>
#include <ShapeImprovementWrapper.hpp>
#include <Mesquite.hpp>

#include <PatchData.hpp>
#include <TerminationCriterion.hpp>
#include <QualityAssessor.hpp>
#include <QualityMetric.hpp>

#include <ConditionNumberQualityMetric.hpp>
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "EdgeLengthRangeQualityMetric.hpp"
#include "Randomize.hpp"

#include <IdealWeightInverseMeanRatio.hpp>
#include <ShapeImprovementWrapper.hpp>
#include <ConditionNumberQualityMetric.hpp>

using namespace Mesquite2;

extern "C" int S3_InverseMeanRatio(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter) {

  MsqError err;
  Smooth3DMesh mon_maillage;

  mon_maillage.initSmooth3DMesh(nb_cells, nb_nodes, nb_nodes_per_cell,
      nodes_number, x, y, z);
  mon_maillage.markNonRelaxedFixed(relax);

  MSQ_CHKERR(err);
  InstructionQueue queue1;
  IdealWeightInverseMeanRatio mean(err);
  MSQ_CHKERR(err);
  LPtoPTemplate obj_func(&mean, 1, err);
  MSQ_CHKERR(err);
  FeasibleNewton pass1(&obj_func);
  //perform optimization globally
  pass1.use_global_patch();

  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit(n_iter);
  pass1.set_outer_termination_criterion(&tc_outer);

  queue1.set_master_quality_improver(&pass1, err);
  MSQ_CHKERR(err);

  queue1.run_instructions(&mon_maillage, err);
  MSQ_CHKERR(err);
  mon_maillage.epilogue(x, y, z, relax); /* On modifie les coordonnées en conséquence */
  return (0);
}
