#include "main.h"
#include "BranchContinuationTask.h"

using namespace std;
using namespace capd;

/**
 * Auxiliary algorithm. For a given parameter 'a' it finds an approximate PO x(a), z(a).
 * 
 * @param seed=(x,0,z,a)
 * @return (newX,0,newZ,a)
 * 
 */ 
capd::IVector findFixedPoint(capd::IVector seed, int threadId){
  auto& pm = RosslerPMPool::dRosslerPMPool->getSolver(threadId).pm;
  pm.getVectorField().setParameter(0,seed[3].mid().leftBound());
  
  capd::DVector u {seed[0].mid().leftBound(),0.,seed[2].mid().leftBound()};
  capd::DMatrix D(3,3);
  for(int i=0;i<10;++i){
    capd::DVector v = pm(u,D);
    D = pm.computeDP(v,D);
    capd::DMatrix A{ {D(1,1)-1,D(1,3)}, {D(3,1),D(3,3)-1} };
    DVector e = matrixAlgorithms::gauss(A,capd::DVector{v[0]-u[0],v[2]-u[2]});
    u[0] -= e[0];
    u[2] -= e[1];
    if(e.euclNorm()<1e-10) break;
  }
  return capd::IVector{u[0],0.,u[2],seed[3]};
}

/**
 * This is an auxiliary data structure that stores some intermediate results
 * of the computation in validation of bifurcation
 */ 
struct BifValidatorData{
  // u0 - candidate for bif point
  // u - box around u
  // lu, ru will be bounds for the curve of fixed points ar the ends of x-range.
  IVector u, lu, ru;
  double deltaZ, deltaA;
  
  BifValidatorData(IVector u, double deltaX, double deltaZ, double deltaA)
    : u(u), deltaZ(deltaZ), deltaA(deltaA)
  {
    this->u[0] += interval(-1,1)*deltaX;
    lu = leftVector(this->u);
    ru = rightVector(this->u);
  }

  /**
   * This is a rigorous routine. It proves the existence of a curve of fixed points
   * for Poincare map parametrized by 
   *   x->(z(x),a(x)) \in (u_z,a_z) + [-dZ,dZ] \times [-dA,dA]
   *  for x \in u_x
   * 
   * If u_x is a single point then it validates the existence of just one fixed point.
   * 
   * @param [in,out] u0 - candidate for fixed point. 
   *    On success, this vector is updated to contain a bound for the fixed point (or curve)
   * @param [out] - on success grad = (z'(x),a'(x))
   * @returns true if curve of fixed point is validated
   */ 
  bool proveFixPoint_fix_x(IVector& u0, IC2PoincareMap& pm, IVector& grad, double dZ, double dA){ 

    // Take a box around u0
    IVector U = u0;
    U[2] += dZ*interval(-1,1);
    U[3] += dA*interval(-1,1);

    // compute f(x,z,a) = (P_x(x,0,z,a)-x,P_z(x,0,z,a)-z)
    C0HOTripleton s0(u0);
    IVector Pu0 = pm(s0);
    IVector f{Pu0[0]-u0[0],Pu0[2]-u0[2]};
    
    // compute derivative of f
    IMatrix D(4,4);
    C1HORect2Set s(U);
    auto v = pm(s,D);
    D = pm.computeDP(v,D);
    IMatrix A{{D(1,3),D(1,4)},{D(3,3)-1,D(3,4)}};

    // compute defect in the interval Newton operator
    IVector N = matrixAlgorithms::gauss(A,f);
    LOGGER(N);

    // additionally compute z'(x) and a'(x)
    grad = -matrixAlgorithms::gauss( A, IVector{D(1,1)-1,D(3,1)} );
    
    // update bound for curve of fixed points (interval Newton opetator)
    u0[2]-=N[0];
    u0[3]-=N[1];
    
    // check inclusion of interval Newton operator and return result
    return subsetInterior(u0[2],U[2]) and subsetInterior(u0[3],U[3]);
  }

  bool proveBifPoint(IC2PoincareMap& c2pm){
    auto& pm = RosslerPMPool::iC2RosslerPMPool->getSolver(0).pm;
    pm.getSolver().setOrder(30);
    IMatrix Dflow(3,3), D(3,3);
    IHessian Hflow(4,4), H(4,4);

    IVector u0 = midVector(u);
    C1HORect2Set s0(IVector{u0[0],u0[1],u0[2]});
    pm.getVectorField().setParameter("a",u0[3]);
    IVector v = pm(s0,Dflow);
    D = pm.computeDP(v,Dflow);
    IVector f{v[0]-u0[0], v[2]-u0[2], (D(1,1)-1)*(D(3,3)-1)- D(1,3)*D(3,1)};

    Dflow = D = IMatrix(4,4);
    IVector U = u0;
    U[0] += interval(-1,1)*3e-11;
    U[2] += interval(-1,1)*1e-14;
    U[3] += interval(-1,1)*6e-14;
    C2Rect2Set s(U);
    v = c2pm(s,Dflow,Hflow);
    c2pm.computeDP(v,Dflow,Hflow,D,H);
    
    interval a = 2*H(0,0,0)*(D(3,3)-1) + (D(1,1)-1)*H(2,0,2) - H(0,0,2)*D(3,1) - 2*D(1,3)*H(2,0,0);
    interval b = 2*H(0,2,2)*(D(1,1)-1) + (D(3,3)-1)*H(2,0,2) - H(0,0,2)*D(1,3) - 2*D(3,1)*H(2,3,3);
    interval c = H(0,1,3)*(D(3,3)-1) + (D(1,1)-1)*H(2,2,3) - H(0,2,3)*D(3,1) - D(1,3)*H(2,0,3);
    IMatrix A{
        {D(1,1)-1,D(1,3),D(1,4)},
        {D(3,1),D(3,3)-1,D(3,4)},
        {a,b,c}
      };
    IVector e = matrixAlgorithms::gauss(A,f);
    LOGGER(e);
    interval eigenvalue = D(3,1) + abs(D(3,3));
    LOGGER(eigenvalue);
    u0[0] -= e[0];
    u0[2] -= e[1];
    u0[3] -= e[2];

    LOGGER(subset(u0,U));
    return subset(u0,U) and abs(eigenvalue)<2e-4;
  }


  /**
   * This is a rigorous routine. 
   * Given a candidate u for a bifurcation point it is shown that
   * 
   * 1. there is a smooth curve of fixed points for Poincare map parametrized by
	   *      x->(z(x),a(x)) \in (u_z,a_z) + [-deltaZ,deltaZ] \times [-deltaA,deltaA]
   *    for x \in u_x + [-deltaX,deltaX] =: [x_*, x^*]
   * 
   * 2.  a'(x_*) and a'(x^*) are of different signs and a''(x)>0 for all x.
   *     Thus a(x) has a uniqe minimum in [x_*, x^*]
   * 
   * @param [in,out] u - approximate bifurcation point. 
   *      On success, this vector is updated to contain a bound for the curve of fixed points
   * @param deltaX, deltaZ, deltaA - sizes of box for validation x->(z(x),a(x)) 
   */ 
  bool proveBifurcation(IC2PoincareMap& pm){
    if(! proveBifPoint(pm) ) return false;

    IVector gradLeft, gradRight, grad;
    // 1. Check the existence of a branch of PO: x->(z(x),a(x))
    if ( !proveFixPoint_fix_x(u,pm,grad,deltaZ,deltaA) ) return false;

    // 2a. Check if a'(x) changes the sign at the ends or x-interval.
    //   Valiadate at the ends of x-range to get tigther bounds for a'
    if( !proveFixPoint_fix_x(lu,pm,gradLeft,7e-8,1e-8) or
        !proveFixPoint_fix_x(ru,pm,gradRight,7e-8,1e-8)) return false;
    // check if the signs are different
    LOGGER(gradLeft);
    LOGGER(gradRight);
    if( !(gradLeft[1]*gradRight[1]<0) ) return false;
    
    // 2b. Check if a''(x)>0
    C2Rect2Set s(u);
    IMatrix D(4,4), Dflow(4,4);
    IHessian H(4,4), Hflow(4,4);
    auto v = pm(s,Dflow,Hflow);
    pm.computeDP(v,Dflow,Hflow,D,H);
    
    // Second derivatives (z''(x),a''(x)) solve this nice 2x2 linear system (implicit equation)
    IMatrix A{{D(1,3),D(1,4)},{D(3,3)-1,D(3,4)}};
    IVector e{H(0,0,0) + H(0,0,2)*grad[0] + (H(0,0,3) + H(0,2,3)*grad[0])*grad[1] + H(0,2,2)*sqr(grad[0]) + H(0,3,3)*sqr(grad[1]),
              H(2,0,0) + H(2,0,2)*grad[0] + (H(2,0,3) + H(2,2,3)*grad[0])*grad[1] + H(2,2,2)*sqr(grad[0]) + H(2,3,3)*sqr(grad[1])};
    IVector hess = -matrixAlgorithms::gauss(A,2*e);
    LOGGER(hess[1]);
    if( !(hess[1]>0) ) return false;

    return true;
  }
};

/**
 * This is an auxiliary algorithm that creates tasks for continuation of a branch of periodic orbits.
 * The branch of fixed point will be parametrized by
 *     a \in [ seed[3].leftBoond(), TrappingRegion::aMax ]
 * We start continuation from the end point of the curve segment obtained in validation of the existence of SN bifurcation.
 * This is
 *     a=seed[3].leftBoond() 
 * @param [in] seed\in R^4 - an approximate fixed point for Poincare map for a parameter value seed[3]
 * @param [in] deltaX, deltaZ sizes of box around an approximate fixed point for the Krawczyk method in BranchContinuationTask
 * 
 * @param [out] list of tasks to be scheduled for execution. These tasks cover entire parameter range [ seed[3].leftBoond(), TrappingRegion::aMax ]
 */ 
void createContinuationTasks(std::vector<BranchContinuationTask*>& tasks, capd::IVector seed, double deltaX, double deltaZ, double deltaA){
  double L = seed[3].leftBound();
  double M = TrappingRegion::aMax;
  for(int i = 0;i<50 and L<M;++i){
    double R = capd::min(M,L+deltaA);
    seed = findFixedPoint(IVector{seed[0],0,seed[2],interval(L,R)},0);
    tasks.push_back( new BranchContinuationTask(seed,deltaX,deltaZ) );
    L = R;
  }
  deltaA *= 10;
  deltaX *= 10;
  deltaZ *= 5;
  while(L<M){
    double R = capd::min(M,L+deltaA);
    seed = findFixedPoint(IVector{seed[0],0,seed[2],interval(L,R)},0);
    tasks.push_back( new BranchContinuationTask(seed,deltaX,deltaZ) );
    L = R;
  }
}

/**
 * This function checks, if segments of periodic orbits validated in different tasks merge into a smooth curve.
 * Here we use the uniqueness property of the Interval Newton Operator.
 * For two subsequent intervals of parameters: 
 *  prev=[a_{i-1},a_{i}] and next=[a_{i},a_{i+1}]
 * we check that a tight bound on the curve at specific parameter 
 *  (x(a_i),z(a_i))
 * belongs to bounds for entire segments prev and next. From the uniqueness we obtain that these two curves must merge at a_i.
 */ 
bool checkContinuity(std::vector<BranchContinuationTask*> branch, capd::IVector u){
  // Sort tasks by increasing parameter 'a'
  std::sort(branch.begin(),branch.end(),[](auto* t1, auto* t2)->bool {return t1->u[3].leftBound() < t2->u[3].leftBound();});
  // start from a bound from parametrization near bifurcation point. This parametrization is x->(x,0,z(x),a(x))
  IVector prev{u[0],u[2]};
  bool result = true;
  for(unsigned int i=0;i<branch.size() and result;++i){
    result = subset(prev[0],branch[i]->u[0]) and subset(prev[1],branch[i]->u[2]);
    prev = branch[i]->Nright;
  }
  return result;
}


/**
 * This is the main algorithm for valiadtion of SN bifurcation and further continuation of two branches of POs.
 * @param [in] bvd - candidate for bifurcation point. In this data structure we store result of subsequent steps of validation procedure.
 * @param [in] deltaX, deltaZ sizes of box around an approximate fixed point for the Krawczyk method in continuation of fixed points. 
*/
bool provePOBranch(BifValidatorData& bvd, IC2PoincareMap& pm, double deltaX, double deltaZ, double deltaA){
  // start from validation of SN bifurcation
  bool result = bvd.proveBifurcation(pm);

  // given curve of fixed points near bifurcation and parametrized by
  // x -> (x,0,z(x),a(x)), x = [x_*,x^*]
  // we create tasks for continuation of two branches of PO
  // they start at 
  //  bvd.lu=(x_*,0,z(x_*),a(x_*)) and bvd.ru=(x^*,0,z(x^*),a(x^*)), respectively 
  std::vector<BranchContinuationTask*> lBranch;
  std::vector<BranchContinuationTask*> rBranch;
  if(result){
    createContinuationTasks(lBranch,bvd.lu,deltaX,deltaZ,deltaA);
    createContinuationTasks(rBranch,bvd.ru,deltaX,deltaZ,deltaA);
  }
  // call general algorithm (include/RosslerPMPool.h") and print result
  int counter = 0;
  auto checkNewtonInclusion = [](BranchContinuationTask* task) { return task->completed  and task->result; };
  result = result and checkCondition(lBranch,checkNewtonInclusion,counter) and checkCondition(rBranch,checkNewtonInclusion,counter);
  // if all tasks return success, we have to check if the segments from each task join into a continuous curve
  result = result and checkContinuity(lBranch,bvd.lu) and checkContinuity(rBranch,bvd.ru);
  LOGGER(counter);
  LOGGER(result);
  return result;
}

int main(int argc, char* argv[])
{
  cout.precision(16); 
  try{
    // Init a pool of objects that compute Poincare map for the Rossler system.
    // The size of the pool=hardware_concurrency is the maximal available number of threads on the machine.
    RosslerPMPool::init(std::thread::hardware_concurrency());
    RosslerPMPool::iRosslerPMPool->setOrder(20);
  
    std::unique_ptr<RosslerPMPool::IC2PMap> PM(RosslerPMPool::IC2PMapFactory("par:b,c;var:x,y,z,a;fun:-(y+z),x+a*y,b+z*(x-c),0;",8,2,capd::ICoordinateSection(4,1),capd::poincare::PlusMinus).createSolver());
    PM->vectorField.setParameter("b",0.2);
    PM->vectorField.setParameter("c",15);

    // Approximate bifurcation points and
    // hand adjusted sizes of set (X,Z,A) on which bifurcation is validated 
    BifValidatorData bvd[] = {
        BifValidatorData(IVector{-24.98615641824929,0,0.005003695953767296,0.2445890212249042}, 5e-4, 3e-7 ,2e-5),
        BifValidatorData(IVector{-27.90101006311546,0,0.004663547021688063,0.3119866509093180}, 2e-4, 3e-7 ,2e-5),
        BifValidatorData(IVector{-29.22211573599117,0,0.004524155280447693,0.3405236989996505}, 1e-4, 5e-8 ,2e-6),
        BifValidatorData(IVector{-29.89377365336397,0,0.004456435009829962,0.3546641965836549}, 4e-5, 5e-8 ,1e-6),
        BifValidatorData(IVector{-30.25817997684925,0,0.004420535075303823,0.3623180183723328}, 1e-5, 6e-9 ,2e-7)
        };

    // Hand adjusted sizes of sets for continuation of periodic orbits for each SN bifurcation
    auto clock = tic("Verification of SN bifurcations and branches of periodic orbits");
      bool result = true;
      result = result and provePOBranch(bvd[0],PM->pm,1e-4,5e-7,1e-6);
      result = result and provePOBranch(bvd[1],PM->pm,1e-5,2e-7,5e-7);
      result = result and provePOBranch(bvd[2],PM->pm,1e-5,2e-7,5e-7);
      result = result and provePOBranch(bvd[3],PM->pm,5e-6,1e-7,5e-7);
      result = result and provePOBranch(bvd[4],PM->pm,1e-6,1e-7,5e-8);
      LOGGER(result);
    tac(clock); // stop clock and print time of computation

  }catch(const exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }

  return 0; 
}
