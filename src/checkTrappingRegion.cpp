///########################################################################
/**
 * @author Daniel Wilczak
 * @date 2025-01-23
 * 
 * This program validates the existence of a trapping region for the Poincare map
 *    P:\Pi\to\Pi,   \Pi = { (x,0,z) : y'<0 }.
 * It is shown that for all parameter values 
 *    a \in [aMin,aMax]
 * P is well defined on the set
 *    W = [xMin,xMax] \times [zMin, zMax]
 * and P(W)\subset W
 * 
 * The constants 
 *    aMin, aMax, xMin, xMax, zMin, zMax 
 * are defined in the file include/TrappingRegion.h as a constexpr static members. 
 * 
 * The program has the following structure:
 * 1) the parameter range [aMin,aMax] is split uniformly into subintervals
 * 2) for each subinterval 'a' the main algorithm validateTrappingRegion is called
 * 3) The main algorithm consists of the following steps (as described in the article)
 * 3a) Cover W =\bigcup_i W_i by boxes and show that P(W_i) is defined (no restrictions on bound of P(W_i) are necessary).
 *     Such computation is performed as a spearate task (instatnce of ComputeImageTask) and executed in a thread pool.
 *     If the existence of P(W_i) cannot be validated, the set W_i is plit into 4 new boxes and we repeat the computation.
 *     If the maximal depth of subdivision (see include/RosslerPMPool.h) is exceeded the the algorithm returns false.
 *     Otherwise, we have validated in final computation that P(W) exists.
 * 3b) We cover boundary of W by a finite number of boxes W_i and we show that P(W_i)\subset W.
 *     Since P is diffeomorphism onto image we conclude that P(W)\subset W.
 *     Again, the computation is heavily parallelized using ComputeImageTask objects and executed in a thread pool.
 */ 
///########################################################################

#include "main.h"
using namespace std;

/**
 * This is a single task to be executed in a thread pool.
 * It computes a bound on Poincare map P(x).
 * On success: 
 * - ComputeImageTask::completed=true 
 * - ComputeImageTask::result = computed bound on P(x)
 * On failure:
 * - ComputeImageTask::completed=false
 * - ComputeImageTask::result = 0
*/ 
struct ComputeImageTask : public capd::threading::Task{
  /// @param _x - an argument for Poincare map 
  ComputeImageTask(const capd::IVector& _x) 
    : x(_x), result(_x), completed(false)
  {}
  
  // async computation of a bound on PoincareMap(x) 
  // this method is executed by a worker in a thread pool 
  void run(unsigned id){    
    try{
      C0HOTripleton set(x); // envelop for box x
      // take a worker from the thread pool and call computation of Poincare map
      result = RosslerPMPool::iRosslerPMPool->getSolver(id).pm(set); 
      completed = true;
    }catch(...){
      result.clear();
    }
  }
  
  // If the image cannot be computed (exception thrown or the bound on the image does not satisfies required properties 
  // then we split the domain x and enqueue new ComputeImage tasks.
  // Subdivision is just the bisection of (x,a) coordinates, so 4 new tasks are created.
  void splitTask(std::vector<ComputeImageTask*>& tasks){
    ComputeImageTask* s1 = new ComputeImageTask(x);
    ComputeImageTask* s2 = new ComputeImageTask(x);
    ComputeImageTask* s3 = new ComputeImageTask(x);
    ComputeImageTask* s4 = new ComputeImageTask(x);
    double xC = (x[0].leftBound()+x[0].rightBound())/2;
    double aC = (x[3].leftBound()+x[3].rightBound())/2;
    s1->x[0].setRightBound(xC); s1->x[3].setRightBound(aC);
    s2->x[0].setLeftBound(xC);  s2->x[3].setRightBound(aC);
    s3->x[0].setRightBound(xC); s3->x[3].setLeftBound(aC);
    s4->x[0].setLeftBound(xC);  s4->x[3].setLeftBound(aC);
    tasks.push_back(s1);
    tasks.push_back(s2);  
    tasks.push_back(s3);
    tasks.push_back(s4);  
  }

  capd::IVector x, result;
  bool completed;
};

int edgesCounter = 0;
int domainCounter = 0;

//########################################################################
/// An auxiliary funciton that creates an initial grid of the trapping region 
/// for further verification if the existence of Poincare map over trapping region.
void createDomainTasks(TrappingRegion& tr, std::vector<ComputeImageTask*>& tasks){
  for(unsigned i=1;i<tr.grid.size();i++)
    tasks.push_back(new ComputeImageTask(tr.gridElement(i-1)));
}

//########################################################################
/// An auxiliary function that creates an initial grid of the boundary of trapping region
/// for further verification if the boundary is mapped into the trapping region.
/// Then, entire trapping region must be mapped into itself.
void createEdgesTasks(TrappingRegion& tr, std::vector<ComputeImageTask*>& tasks){
  for(unsigned i=1;i<tr.grid.size();i++){
    IVector u = tr.gridElement(i-1); // a set of the form [x]\times [0] \times [z]
    IVector u1=u, u2=u;
    u1[2] = u1[2].rightBound(); // u1 = [x]\times [0] \times [zMaxz]
    u2[2] = u2[2].leftBound();  // u1 = [x]\times [0] \times [zMin]
    tasks.push_back(new ComputeImageTask(u1));
    tasks.push_back(new ComputeImageTask(u2));
  }
  // finally, add left and right edges 
  // [min(x)] \times [0] \times [z]
  // [max(x)] \times [0] \times [z]
  tasks.push_back(new ComputeImageTask(tr.getLeftEdge()));
  tasks.push_back(new ComputeImageTask(tr.getRightEdge()));
}

//########################################################################
/**
 * This is the main algorithm for validation that the rectangle 
 * described by data structure TrappingRegion is indeed a trapping region for Poincare map.
 * This is a two-step verification.
 * 1. Check if the Poincare map indeed exists for every element in the grid.
 * 2. Check if the boundary of the trapping region (rectangle) is mapped into the itself.
 */  
bool validateTrappingRegion(interval a)
{
  // predicates for both steps
  auto doesPoincareMapExist = [](ComputeImageTask* task) { return task->completed; };
  auto isEdgeMappedIntoTrappingRegion = [](ComputeImageTask* task) { 
      IVector u = task->result;
      return task->completed  and
             u[0].leftBound()  >= TrappingRegion::xMin and 
             u[0].rightBound() <= TrappingRegion::xMax and
             u[2].leftBound()  >= TrappingRegion::zMin and 
             u[2].rightBound() <= TrappingRegion::zMax;    
  };

  // create initial grids of the trapping region and its boundary
  TrappingRegion tr(a);
  std::vector<ComputeImageTask*> domainTasks, edgesTasks;
  createDomainTasks(tr,domainTasks);
  createEdgesTasks(tr,edgesTasks);

  // call abstract algoritm from checkCondition from include/RosslerPMPool.h
  // that checks if 'condition' is satisfied for all tasks
  bool result = true;
  result = result and checkCondition(domainTasks,doesPoincareMapExist,edgesCounter);
  result = result and checkCondition(edgesTasks,isEdgeMappedIntoTrappingRegion,domainCounter);
  for(ComputeImageTask* task : domainTasks) delete task;
  for(ComputeImageTask* task : edgesTasks) delete task;

  return result;	
}

//########################################################################

int main(int argc, char* argv[])
{
  cout.precision(16); 
  try{
    // Init a pool of objects that compute Poincare map for the Rossler system.
    // The size of the pool=hardware_concurrency is the maximal available number of threads on the machine.
    RosslerPMPool::init(std::thread::hardware_concurrency());
    RosslerPMPool::iRosslerPMPool->setOrder(15);
    RosslerPMPool::iRosslerPMPool->setTolerance(1e-7,1e-7);

    interval da = interval(TrappingRegion::aMax) - interval(TrappingRegion::aMin);
    bool result = true;

    auto clock = tic("Validation of the existence of a trapping region");
      // Split parameter range [aMin,aMax] into N subintervals and check the existence of trapping region for each subinterval. 
      constexpr int N = 2048;
      for(int i=0;i<N and result;i++){
        interval a = TrappingRegion::aMin + da*interval(i,i+1)/N;
        LOGGER(i);
        LOGGER(a);
        // here we call a routine that validates P(W)\subset W for a subinterval of parameters
        result = result and validateTrappingRegion(a);
        std::cout << "#######################\n\n";
      }
      LOGGER(domainCounter);
      LOGGER(edgesCounter);
      LOGGER(result);
    tac(clock); // stop clock and print time of computation
  }catch(const exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }
  return 0;
}
