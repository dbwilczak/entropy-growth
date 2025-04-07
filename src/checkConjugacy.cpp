///########################################################################
/**
 * @author Daniel Wilczak
 * @date 2025-01-23
 * 
 * This program constructs and validates semiconjugacy of the Poincare map
 *    P:\Pi\to\Pi,   \Pi = { (x,0,z) : y'<0 }
 * and restricted to the trapping region W=[xMin,xMax]\tiems [zMin,zMax] to a shift dynamics.
 *
 * The constants of the trapping region and the parameter range
 *    aMin, aMax, xMin, xMax, zMin, zMax 
 * are defined in the file include/TrappingRegion.h as a constexpr static members. 
 * 
 * 
 * The program has the following structure:
 * 1) The parameter range [aMin,aMax] is initialy split uniformly into subintervals.
 * 2) For each subinterval 'a' the main algorithm checkConjugacy is called.
 * 3) The main algorithm consists of the following steps (as described in the article)
 * 3a) nonrigorous routine: detect approximate extrema of the function x -> P_x(x,0,z), for z=zMin, zMax and zMid = (zMin+zMax)/2
 * 3b) nonrigorous routine: based on the above extrema predict location of h-sets for verification of covering relations
 * 3c) start CheckConjugacyTask which compute bounds on Poincare map and validate transition matrix for shift dynamics.
 * 4) Finally, subintervals with the same transition matrix are merged and the final transition matrices are printed.  
 */ 
///########################################################################

#include <sstream>
#include "main.h"
#include "ConjugacyData.h"
using namespace std;

/**
 * This is a single task to be executed in a thread pool.
 * 
 * The virtual method 'run' executed by a thread pool 
 * 1) calls subroutine 
 *      this->constructConjugacy(); 
 *    to contruct h-sets and a guess for transition matrix
 * 2) calls subroutine
 *      this->checkConjugacy();
 *    This function comutes bounds on P(e), e\in E, 
 *    where E is finite set of edges of constructed h-sets.
 *    Finally condition for covering relation is checked and 1,0,-1 
 *    are written to conjugacy matrix in a proper place.
 * 3) calls subroutine
        this->removeEmptyRows();
      It removes empty (zero) rows and corresponding columns from the transition matrix.
      Such a row may appear if some predicted covering relation could not be validated.
*/

/// Base class ConjugacyData contains implementation of nonrigorous algorithms for prediction of symbolic dynamics,
/// as well as data fro symbolic dynamics, such as extrema and edges of h-sets.
struct CheckConjugacyTask : public capd::threading::Task, public ConjugacyData 
{
  bool result = false;
  bool completed = false;
  
  // Store subinterval a and some threshold parameter m>0 used to construct h-sets
  CheckConjugacyTask(interval a, double m) : ConjugacyData(a,m) {}

  /// Auxiliary routine - computes Poincare map P_x(u)
  interval eval(const IVector& u){
    C0HOTripleton s(u);
    return RosslerPMPool::iRosslerPMPool->getSolver(id).pm(s)[0];
  }

  /// Auxiliary routine - compute P_x(x,0,[zMin,zMax])
  interval eval(std::vector<IVector>& grid, int i, int depth, int dir, std::map<double,interval>& computed){
    double _x = grid[0][0].leftBound();
    auto it = computed.find(_x);
    if(it!=computed.end()) return it->second;
    
    std::vector<IVector> tmp;
    interval rX = 0.;
    while(depth and grid.size()){
      for(auto& u : grid){
        interval x = eval(u);
        bool result = true;
        for(unsigned j=0;j<symbols.size() and result;++j){
          if(conjugacy(i+1,j+1)==0) continue;
          if(conjugacy(i+1,j+1)==1 and dir==1 and x < symbols[j].first) continue;
          if(conjugacy(i+1,j+1)==1 and dir==-1 and x > symbols[j].second) continue;
          if(conjugacy(i+1,j+1)==-1 and dir==-1 and x < symbols[j].first) continue;
          if(conjugacy(i+1,j+1)==-1 and dir==1 and x > symbols[j].second) continue;
          result = false;
        }
        if(!result){ // split on failure
          tmp.push_back(u);
        } else {
          if(rX==0.)
            rX = x;
          else
            rX = intervalHull(rX,x);
        }
      }
      grid.clear();
      for(auto& u : tmp){
        double _z = (u[2].leftBound()+u[2].rightBound())/2;
        double _a = (u[3].leftBound()+u[3].rightBound())/2;
        IVector u1 = u;
        IVector u2 = u;
        IVector u3 = u;
        IVector u4 = u;
        u1[2].setLeftBound(_z); u1[3].setLeftBound(_a);
        u2[2].setLeftBound(_z); u2[3].setRightBound(_a);
        u3[2].setRightBound(_z); u3[3].setLeftBound(_a);
        u4[2].setRightBound(_z); u4[3].setRightBound(_a);
        grid.push_back(u1);
        grid.push_back(u2);
        grid.push_back(u3);
        grid.push_back(u4);
      }
      depth--;
    }
    computed[_x] = rX;
    return rX;
  }

  /// Compute bounds on edges (using eval algorithm) and check inequalities for covering relation
  void checkConjugacy(){
    std::map<double,interval> computed;
    this->result = true;
    int depth = 20, N = 1;
    interval d = (interval(TrappingRegion::zMax)-interval(TrappingRegion::zMin))/N;

    for(unsigned i=0;i<this->symbols.size() and this->result;++i){
      bool fresult = true;
      bool sresult = true;
      std::vector<IVector> fgrid, sgrid;
      for(int k=0;k<N;++k){
        fgrid.push_back( IVector{symbols[i].first,0.,TrappingRegion::zMin + interval(k,k+1)*d,a} );
        sgrid.push_back( IVector{symbols[i].second,0.,TrappingRegion::zMin + interval(k,k+1)*d,a});
      }
      // take two edges of a set and compute bounds using eval algorithm
      interval x1 = eval(fgrid,i,depth,1,computed);
      interval x2 = eval(sgrid,i,depth,-1,computed);
      // check if inequalities for covering relation hold true
      for(unsigned j=0;j<symbols.size() and result;++j){
        if(conjugacy(i+1,j+1)==1)
          result = result and x1 < symbols[j].first and x2>symbols[j].second;
        if(conjugacy(i+1,j+1)==-1)
          result = result and x2 < symbols[j].first and x1>symbols[j].second;
        // if not, set zero to transition matrix
        if(!result){
          conjugacy(i+1,j+1) = 0;
          result = true;
        }
      }
    }
  }
  
  /// Auxiliary algorithm - remove empty rows/columns from transition matrix
  void removeEmptyRows(){
    for(int i=1;i<=conjugacy.numberOfRows();++i){
      int sum = 0;
      for(int m=1;m<=conjugacy.numberOfRows();++m)
        sum += conjugacy(i,m);
      if(sum==0){
        ZMatrix m(conjugacy.numberOfRows()-1,conjugacy.numberOfRows()-1);
        int j,r;
        for(j=1,r=1;j<=conjugacy.numberOfRows();++j){
          if(j==i) continue;
          int k,s;
          for(k=1,s=1;k<=conjugacy.numberOfRows();++k){
            if(k==i) continue;
            m(r,s) = conjugacy(j,k);
            ++s;
          }
          ++r;
        }
        swap(conjugacy,m);
      }
    }
  }

  // Async construction and verification of semiconjugacy to shift dynamics.
  // This method is executed by a worker in a thread pool.
  void run(unsigned id){    
    this->id = id; // save threadId
    try{
      this->ConjugacyData::constructConjugacy();
      this->checkConjugacy();
      this->removeEmptyRows();
    }catch(std::exception& e){
      result = false;
    }
    completed = true;
  }  
};

// ConjugacyResult is a list of pairs (interval,transition matrix)
typedef std::vector< std::pair<interval,ZMatrix> > ConjugacyResult;

//########################################################################

/**
 * This is the main algorithm
 * @param [in] s - size of subinterval of the parameter range [aMin,aMax]
 * @param [in] margin>0 - a control parameter used to construct h-sets
 * @param [out] computed list of pairs (interval,transition matrix)
 * @returns true if computation are completed
 *          false if come exception is thrown during the computation
 */ 
bool checkConjugacy(double s, double margin, ConjugacyResult& r){
  // create a list of task, that is subdivide [aMin,aMax]
  std::vector<CheckConjugacyTask*> tasks;
  double L = TrappingRegion::aMin;
  double M = TrappingRegion::aMax;
  do{
    double R = capd::min(L+s,M); 
    interval a(L,R);
    tasks.push_back(new CheckConjugacyTask(a,margin));
    L = R;
  }while(L<M);

  // Execute and join tasks
  for(CheckConjugacyTask* task : tasks) RosslerPMPool::threadPool->process(task);
  for(CheckConjugacyTask* task : tasks) {
    task->join();
  }
  
  // Check, if all tasks completed the computation and returned result
  // Collect data from tasks
  bool result = true;
  for(CheckConjugacyTask* task : tasks) {
    result = result and task->result;
    if(result) r.push_back(std::pair(task->a,task->conjugacy));
    delete task;
  }
  return result;
}

//########################################################################
/// Auxiliary function that prints result of computation in the LaTeX format.
void mergeAndPrintResult(ConjugacyResult& cr){
  ConjugacyResult cr_merged;

  interval a = cr[0].first;
  ZMatrix M = cr[0].second;
  for(unsigned i=1;i<cr.size();++i){
    if( cr[i].second.dimension() == M.dimension() and  cr[i].second==M ) 
      a = intervalHull(a,cr[i].first);
    else {
      cr_merged.push_back(std::pair(a,M));
      a = cr[i].first;
      M = cr[i].second;
    }
  }
  cr_merged.push_back(std::pair(a,M));
  
  // print conjugacy matrices (further computation of max eigenvalue in external tool)
  for(auto& r : cr_merged){
    auto [a,M] = r;
    for(auto& i : M) i = capd::abs(i);
    std::cout << a << "\n" << M << "\n\n";
  }
}

//########################################################################

int main(int argc, char* argv[]){
  std::cout.precision(17);
  try{
    // Init a pool of objects that compute Poincare map for the Rossler system.
    // The size of the pool=hardware_concurrency is the maximal available number of threads on the machine.
    RosslerPMPool::init(std::thread::hardware_concurrency());
    RosslerPMPool::iRosslerPMPool->setOrder(20);
    RosslerPMPool::iRosslerPMPool->setTolerance(1e-21,1e-21);

    auto clock = tic("Validation of semiconjugacy to shift dynamics");
      // Call the main algorithm. Split parameter range into subintervals of width 1./1<<20 = 2^{-20} \approx 1e-6.
      // tol=1e-3 is a control parameter used to predict h-sets for symbolic dynamics.
      ConjugacyResult cr;
      bool result = checkConjugacy(1e-6,1e-3,cr);

      // merge subintervals of parameter range with the same transition matrix and print result
      if(result)
        mergeAndPrintResult(cr);
      LOGGER(result);
    tac(clock); // stop clock and print time of computation
  }catch(const exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }
  return 0;
}
