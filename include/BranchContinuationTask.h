#ifndef ROSSLER_TEMPLAES_BRANCH_CONTINUATION_TASK_H
#define ROSSLER_TEMPLAES_BRANCH_CONTINUATION_TASK_H

#include "main.h"
#include "capd/dynset/C0HOSet.hpp"
#include "capd/dynset/C0TripletonSet.hpp"

capd::IVector findFixedPoint(capd::IVector seed, int threadId);

/**
 * This is a single task to be executed in a thread pool.
 * It validates the existence of a segment of fixed points of the Poincare map
 *    a->u(a), P(u(a))=u(a)
 * parametrized by a =[a_*,a^*] = _u[3].
 * The algorithm applies interval Newton method to the function 
 *    f(u) = P(u)-u
 * On success: completed=true and the bound on this curve is stored in member N.
 */ 
struct BranchContinuationTask : public capd::threading::Task{
  BranchContinuationTask(const capd::IVector& _u, double dX, double dZ) 
    : u(_u), Pu0(_u), N(_u), Nright(_u), deltaX(dX), deltaZ(dZ)
  {
    u[0] += deltaX*interval(-1,1);
    u[2] += deltaZ*interval(-1,1); 
  }

  IVector computeNewtonOperator(capd::IPoincareMap& pm, IVector X){
      // compute P(x0)
      capd::IVector u0 = midVector(X);
      u0[3] = X[3];
      C0HOTripleton s0(u0);
      Pu0 = pm(s0);

      // compute DP(x)
      C1HORect2Set s(X);
      capd::IMatrix D(4,4);
      capd::IVector v = pm(s,D);
      D = pm.computeDP(v,D);

      // compute interval Newton operator
      capd::IVector f{Pu0[0]-u0[0],Pu0[2]-u0[2]};
      capd::IMatrix A{{D(1,1)-1,D(1,3)},{D(3,1),D(3,3)-1}};
      return capd::IVector{u0[0],u0[2]} - capd::matrixAlgorithms::gauss(A,f);
  }

  // async validation of the existence of a branch of POs parametrized by 'a'.
  void run(unsigned id){    
    try{
      auto& pm = RosslerPMPool::iRosslerPMPool->getSolver(id).pm;
      this->N = computeNewtonOperator(pm,u);
      result = subsetInterior(N[0],u[0]) and subsetInterior(N[1],u[2]);
      // If validated, compute addidtionally tight bound on the fixed point at the end of interval 'a'
      // This is needed to validate that subsequent segments merge into a continuous curve.
      if(result){
        IVector X = u;
        X[3] = u[3].rightBound();
        X = findFixedPoint(X,id);
        // hardcoded size of set for interval Newton operator at the end of interval of parameters
        // this is needed to check, that the curve of fixed points 
        // validated in different subintervals of parameter glue to a continuous curve
        X[0] += 1e-7*interval(-1,1);
        X[2] += 1e-7*interval(-1,1);
        this->Nright = computeNewtonOperator(pm,X);
        result = subsetInterior(Nright[0],X[0]) and subsetInterior(Nright[1],X[2]);
        this->Nright = intersection(this->Nright,this->N);
      }
      completed = true;     
    }catch(...){
      N.clear();
      Nright.clear();
      Pu0.clear();
    }

    // on failure, prepare data for split of the range of parameters
    if(! (completed and result) ){
      double aC = (u[3].leftBound()+u[3].rightBound())/2;
      uL = u;
      uR = u;
      uL[3].setRightBound(aC);
      uR[3].setLeftBound(aC);
      uL = findFixedPoint(uL,id);
      uR = findFixedPoint(uR,id);
    }
  }
  
  // If image cannot be computed or the bound on the image does not satisfies required properties 
  // then we split the domain x and enqueue new ComputeImage tasks.
  // Subdivision is just the bisection of the parameter range 'a', so 2 new tasks are created.
  void splitTask(std::vector<BranchContinuationTask*>& tasks){
    tasks.push_back(new BranchContinuationTask(uL,deltaX,deltaZ));
    tasks.push_back(new BranchContinuationTask(uR,deltaX,deltaZ));  
  }

  capd::IVector u, Pu0, N, Nright;
  capd::IVector uL, uR;
  double deltaX, deltaZ;
  bool completed = false, result = false;
};

#endif
