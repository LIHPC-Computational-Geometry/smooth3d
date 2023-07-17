/*
 * ConditionNumber.cpp
 *
 *  Created on: 21 févr. 2013
 *      Author: weilljc
 *
 *      @Id: @
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

#include <IdealWeightMeanRatio.hpp>
#include <ShapeImprovementWrapper.hpp>
#include <ConditionNumberQualityMetric.hpp>

using namespace Mesquite2;

extern "C" int S3_ConditionNumber(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter) {

  MsqError err;
  Smooth3DMesh mon_maillage;

  mon_maillage.initSmooth3DMesh(nb_cells, nb_nodes, nb_nodes_per_cell,
      nodes_number, x, y, z);
  mon_maillage.markNonRelaxedFixed(relax);

  MSQ_CHKERR(err);

  // creates a mean ratio quality metric ...
  ConditionNumberQualityMetric shape_metric;
  LPtoPTemplate obj_func2(&shape_metric, 2, err);
  ConjugateGradient pass2(&obj_func2, err);

  TerminationCriterion sc3;
  sc3.add_iteration_limit(n_iter);

  pass2.set_inner_termination_criterion(&sc3);

  // creates an intruction queue

  InstructionQueue queue1;
  queue1.set_master_quality_improver(&pass2, err);
  MSQ_CHKERR(err);

  queue1.run_instructions(&mon_maillage, err);
  MSQ_CHKERR(err);

  mon_maillage.epilogue(x, y, z, relax); /* On modifie les coordonnées en conséquence */

  return 0;
}
