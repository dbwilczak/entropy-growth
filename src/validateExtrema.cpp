///########################################################################
/**
 * @author Daniel Wilczak
 * @date 2025-01-23
 * 
 * This program computed exact number of extrema of the function
 *    f(x) := P_x(x,0,z), z\in[zMin,zMax], x\in[xMin,xMax]
 * where P is a Poincare map:
 *    P:\Pi\to\Pi,   \Pi = { (x,0,z) : y'<0 }
 * and for selected parameter values a\in[aMin,aMax].
 *
 * The constants 
 *    aMin, aMax, xMin, xMax, zMin, zMax 
 * are defined in the file include/TrappingRegion.h as a constexpr static members. 
 * 
 * The program has the following structure:
 * 1) In the main function we call a general algorithm 'checkTemplate'
 *    for selected parameter values a\in[aMin,aMax].
 * 2) For each selected parameter 'a' the main algorithm checkTemplate
 *    splits the domain \bigcup_i W_i = X_i \times [zMin,zMax] 
 * 3a) For each set W_i a task computing bound on the derivative f'(W_i) and/or f"(W_i) is created.
 * 3b) If derivatives cannot be comuted (exception thrown) or both f' and f" can be zero, the domain W_i is split. 
 * 4) Finally, the data from completed tasks is collected 
 *    and the exact number of extrema is computed by checking signs of f' and f".
 * 
 */ 
///########################################################################

#include "main.h"
#include "ConjugacyData.h"
using namespace std;

/**
 * This is a simple data struture that stores an argument u=(x,0,z) 
 * and bounds on f'(x) and/or f"(x), where f(x) = P_x(u)
 */ 
struct DerivativeData {
  DerivativeData(capd::IVector u) : u(u) {}
  
  capd::IVector u;
  interval dx = 0.;
  interval d2x = 0.;
  capd::interval Px = 0.;
  
  void print(){
    std::cout.precision(17);
    std::cout << "$"<< u[0] << "$ & $" << Px << "$ & $" << dx << "$ & ";
    if(isSingular(d2x)) std::cout << "$-$";
    else  std::cout << "$" << d2x << "$";
    std::cout << "\\\\ \\hline\n";
  }
  void merge(const DerivativeData& data){
    this->u = intervalHull(this->u,data.u);
    this->Px = intervalHull(this->Px,data.Px);
    this->dx = intervalHull(this->dx,data.dx);
    this->d2x = intervalHull(this->d2x,data.d2x);
  }
};

/**
 * This is a single task to be executed in a thread pool.
 * It computes a bound on the first and/or second derivative of 
 *   f(x) := P(x,0,z) for all z\in[zMin,zMax].
 * On success: 
 * - ComputeDerivativeTask::completed=true 
 * - ComputeDerivativeTask::Px = computed bound on f(x)
 * - ComputeDerivativeTask::dx = computed bound on f'(x)
 * - if deg==2 then ComputeDerivativeTask::d2x = computed bound on f"(x)
 * On failure:
 * - ComputeImageTask::completed=false
 * - ComputeImageTask::result = 0
*/ 
struct ComputeDerivativeTask : public DerivativeData, public capd::threading::Task{
  ComputeDerivativeTask(const capd::IVector& _u, int deg) 
    : DerivativeData(_u), deg(deg)
  {}
  
  // async computation of a bound f(x), f'(x) and/or f"(x) 
  // this method is executed by a worker in a thread pool 
  void run(unsigned id){
    try{
      auto& pm = RosslerPMPool::iC2RosslerPMPool->getSolver(id).pm;

      // split Z=[zMin,zMax] into K subintervals
      constexpr int K = 40;
      std::vector<interval> zGrid, PxGrid, dxGrid, d2xGrid;
      for(int i=0;i<K;++i)
        zGrid.push_back(u[2].leftBound() + interval(i,i+1)*(u[2].right()-u[2].left())/K);
      
      // for each element in zGrid 
      for(auto z:zGrid){
        capd::IVector v({u[0],0.,z});
        if(deg==1){ // enclose derivative
          capd::C1Rect2Set s(v);    
          capd::IMatrix DP(3,3);
          capd::IVector Pu = pm(s,DP);
          PxGrid.push_back(Pu[0]);
          dxGrid.push_back(pm.computeDP(Pu,DP)(1,1));
          d2xGrid.push_back(0.0);
        }else{ // enclode 1st and 2nd derivative
          capd::C2Rect2Set s(v);
          capd::IMatrix  Df(3,3), DP(3,3);
          capd::IHessian Hf(3,3), HP(3,3);
          capd::IVector Pu = pm(s,Df,Hf);
          pm.computeDP(Pu,Df,Hf,DP,HP);
          PxGrid.push_back(Pu[0]);
          dxGrid.push_back(DP(1,1));
          d2xGrid.push_back(HP(0,0,0));
        }
      }
      // merge results from K subintervals
      Px = *PxGrid.begin();  for(auto i : PxGrid) Px = intervalHull(Px,i);
      dx = *dxGrid.begin();  for(auto i : dxGrid) dx = intervalHull(dx,i);
      d2x = *d2xGrid.begin();  for(auto i : d2xGrid) d2x = intervalHull(d2x,i);
      completed = true;
    }catch(...){
      Px=dx=d2x=0.;
      completed = false;
    }
  }
  
  // If derivativces cannot be computed or the bound on the image does not satisfies required properties 
  // then we split the domain x and enqueue two new ComputeDerivativeTasks.
  void splitTask(std::vector<ComputeDerivativeTask*>& tasks){
    ComputeDerivativeTask* s1 = new ComputeDerivativeTask(u,deg);
    ComputeDerivativeTask* s2 = new ComputeDerivativeTask(u,deg);
    double x = (u[0].leftBound()+u[0].rightBound())/2;
    s1->u[0].setRightBound(x);
    s2->u[0].setLeftBound(x);
    tasks.push_back(s1);
    tasks.push_back(s2);
  }

  bool completed = false;
  int deg;
};

/// This routine splits domain [xMin,xMax] into subintervals X_i
/// For each X_i a task that computes f'(x) and/or f"(x) is created
void createTasks(interval a, std::vector<ComputeDerivativeTask*>& tasks){
  // first compute approximate extrema
  std::vector<double> extrema;
  Extrema::findExtrema(a.mid().leftBound(),TrappingRegion::zMid,extrema,0);
  std::sort(extrema.begin(),extrema.end());
  std::cout << "Parameter value a=" << a << std::endl;
  std::cout << "Approximate extrema: ";
  for(double e: extrema)  std::cout << e << " ";
  std::cout << std::endl;
  
  // z = z-range of the trapping region
  interval z = interval(TrappingRegion::zMin,TrappingRegion::zMax);

  // split x-range
  // between extrema we request 1st derivative only,
  // close to approximate extrema the 2nd derivative is requested
  double L = TrappingRegion::xMin;
  for(double e: extrema){   
    double delta = (e-L)/256;
    double maxR = e-delta;    
    // case of 1st derivative
    while(L<maxR){
      double R = capd::min(maxR,L+delta);
      interval x = interval(L,R);    
      tasks.push_back( new ComputeDerivativeTask(IVector({x,0.,z}),1) );      
      L = R;
    }
    // case of 2nd derivative
    double R = e+delta;
    interval x = interval(L,R);    
    tasks.push_back( new ComputeDerivativeTask(IVector({x,0.,z}),2) );      
    L = R;    
  }
  // split final subinterval of x-range between last extremum and xMax
  double delta = (TrappingRegion::xMax-L)/128;
  while(L<TrappingRegion::xMax){
    double R = capd::min(TrappingRegion::xMax,L+delta);
    interval x = interval(L,R);    
    tasks.push_back( new ComputeDerivativeTask(IVector({x,0.,z}),1) );      
    L = R;
  }
}

/// This is the main function that computes exact numbers of extrema 
/// of the function
///    f(x) = P_x(x,0,Z)
///  for a given parameter 'a'. 
void computeExtrema(interval a){
  RosslerPMPool::dRosslerPMPool->setParameter("a",a.mid().leftBound()); 
  RosslerPMPool::iC2RosslerPMPool->setParameter("a",a); 

  // create tasks that compute derivatives
  std::vector<ComputeDerivativeTask*> tasks;
  createTasks(a,tasks);
  
  // we accept result of a task if either f' or f'' is non-vanishing
  // otherwise the task will be split
  auto condition = [](ComputeDerivativeTask* t) -> bool{
    return t->completed and ( !isSingular(t->dx) or !isSingular(t->d2x) );
  };
  // run all tasks (see include/RosslerPMPool.h
  int counter = 0;
  bool result = checkCondition(tasks,condition,counter);
  
  // Sort tasks by increasing x-coordinate of the argument
  std::sort(tasks.begin(),tasks.end(),[](auto* t1, auto* t2)->bool {return t1->u[0].leftBound() < t2->u[0].leftBound();});

  // Merge results and print output in LaTeX format
  bool first = true;
  DerivativeData derData = *tasks.front();
  DerivativeData* prev = nullptr;
  std::vector<DerivativeData> vdd;
  
  for(DerivativeData* task : tasks){
    if(first){
      derData = *task;
      first = false;
      prev = task;
    }
    else{
      if(
        (isSingular(task->dx) and isSingular(derData.dx)) or 
        (!isSingular(task->dx) and !isSingular(derData.dx)) or
        (!isSingular(task->d2x) and !isSingular(derData.d2x)) or
        (task->u[0].leftBound() == prev->u[0].leftBound())
      ){
        derData.merge(*task);
      } else {
        vdd.push_back(derData);
        // reset to a new branch
        derData = *task;
      }
    }
    prev = task;
  }
  vdd.push_back(derData);

  for(ComputeDerivativeTask* task : tasks) delete task;

  // Compute relevant extrema
  std::cout << "Bounds on 1st and/or 2nd derivatives: (x, f(x), f'(x), f\"(x))\n";
  interval eMax = TrappingRegion::xMin, eMin = TrappingRegion::xMax;
  for(DerivativeData& dd : vdd) {
    dd.print(); // print LaTeX table with bounds on derivatives
    if(isSingular(dd.dx)){
      eMax = capd::max(eMax,dd.Px);
      eMin = capd::min(eMin,dd.Px);
    }
  }
  std::vector<DerivativeData> relevantExtrema;
  std::vector<DerivativeData> possibleRelevantExtrema;
  std::copy_if(vdd.begin(),vdd.end(),std::back_inserter(possibleRelevantExtrema),[eMin,eMax](auto& dd) {return isSingular(dd.dx) and dd.u[0]<=eMax.rightBound() and dd.u[0]>=eMin.leftBound();});
  std::copy_if(vdd.begin(),vdd.end(),std::back_inserter(relevantExtrema),[eMin,eMax](auto& dd) {return isSingular(dd.dx) and dd.u[0]<=eMax.leftBound() and dd.u[0]>=eMin.rightBound();});
  result = result and (relevantExtrema.size()==possibleRelevantExtrema.size());
  LOGGER(result);
  if(result){
    std::cout << "Relevant extrema (" << relevantExtrema.size() << "):\n";
    for(DerivativeData& re : relevantExtrema) std::cout << re.u[0] << " ";
    std::cout << std::endl;
  }
  std::cout << "\n###############################\n\n";
}

//########################################################################

int main(int argc, char* argv[]){
 
  cout.precision(17); 
  try{
    // Init a pool of objects that compute Poincare map for the Rossler system.
    // The size of the pool=hardware_concurrency is the maximal available number of threads on the machine.
    RosslerPMPool::init(std::thread::hardware_concurrency());

    // compute exact number of extrema of 
    //    x->P_x(x,0,Z)
    // for selected parameter values, as described in the article
    auto clock = tic("Computation of number of extrema");
      for(int a : {1200, 2000, 2600, 3000, 3200, 3330, 3450, 3500, 3560, 3600, 3620, 3659} ) // we want interval bounds of 0.12, 0.2, etc 
        computeExtrema(interval(a)/10000);
    tac(clock); // stop clock and print time of computation
	
  }catch(const exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }
  return 0;
}
