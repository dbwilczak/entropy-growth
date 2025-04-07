#ifndef ROSSLER_TEMPLAES_ROSSLER_PM_POOL_H
#define ROSSLER_TEMPLAES_ROSSLER_PM_POOL_H

#include <memory>
#include <thread>

#include "capd/capdlib.h"
#include "capd/threading/SolverPool.h"
#include "capd/threading/ThreadPool.h"
#include "capd/threading/SolverFactory.h"

/**
 * This is a class that creates and stores a static pool of instances of PoincareMap objects for the Rossler system.
 * This pool is then used in concurrent algorithms that check some properties of the Rossler system.
*/
struct RosslerPMPool{
  typedef capd::threading::PMap<capd::IPoincareMap,capd::ICoordinateSection> IPMap;
  typedef capd::threading::PMap<capd::DC2PoincareMap,capd::DCoordinateSection> DPMap;
  typedef capd::threading::PMap<capd::IC2PoincareMap,capd::ICoordinateSection> IC2PMap;

  typedef capd::threading::SolverPool<IPMap> IPMPool;
  typedef capd::threading::SolverPool<DPMap> DPMPool;
  typedef capd::threading::SolverPool<IC2PMap> IC2PMPool;
  
  typedef capd::threading::PMapFactory<capd::IPoincareMap,capd::ICoordinateSection> IPMapFactory;
  typedef capd::threading::PMapFactory<capd::DC2PoincareMap,capd::DCoordinateSection> DPMapFactory;
  typedef capd::threading::PMapFactory<capd::IC2PoincareMap,capd::ICoordinateSection> IC2PMapFactory;


  // the main thread pool of working threads that execute async tasks
  static std::unique_ptr<capd::threading::ThreadPool> threadPool;
  
  // instances of Poincare maps for rigorous computation
  static std::unique_ptr<IPMPool> iRosslerPMPool;
  // instances of Poincare maps for non-rigorous computation
  static std::unique_ptr<DPMPool> dRosslerPMPool;
  // instances of Poincare maps for rigorous computation of 2nd order variational equations
  static std::unique_ptr<IC2PMPool> iC2RosslerPMPool;
  
  static void init(int threadsNo){
    threadPool.reset(new capd::threading::ThreadPool(threadsNo));

    // vectorField, order, degree, section, CrossingDirection
    DPMapFactory dFactory("par:a,b,c;var:x,y,z;fun:-(y+z),x+a*y,b+z*(x-c);",20,1,capd::DCoordinateSection(3,1),capd::poincare::PlusMinus);
    IPMapFactory iFactory("par:b,c;var:x,y,z,a;fun:-(y+z),x+a*y,b+z*(x-c),0;",15,1,capd::ICoordinateSection(4,1),capd::poincare::PlusMinus);
    IC2PMapFactory iC2Factory("par:a,b,c;var:x,y,z;fun:-(y+z),x+a*y,b+z*(x-c);",15,2,capd::ICoordinateSection(3,1),capd::poincare::PlusMinus);

    dRosslerPMPool.reset(new DPMPool(threadsNo,dFactory));
    iRosslerPMPool.reset(new IPMPool(threadsNo,iFactory));
    iC2RosslerPMPool.reset(new IC2PMPool(threadsNo,iC2Factory));

    dRosslerPMPool->setParameter("b",0.2);
    iRosslerPMPool->setParameter("b",interval(2)/10);
    iC2RosslerPMPool->setParameter("b",interval(2)/10);

    dRosslerPMPool->setParameter("c",15);
    iRosslerPMPool->setParameter("c",15);
    iC2RosslerPMPool->setParameter("c",15);    
  }
};

//########################################################################
/// Definitions of static members
std::unique_ptr<RosslerPMPool::IPMPool> RosslerPMPool::iRosslerPMPool; 
std::unique_ptr<RosslerPMPool::DPMPool> RosslerPMPool::dRosslerPMPool;
std::unique_ptr<RosslerPMPool::IC2PMPool> RosslerPMPool::iC2RosslerPMPool; 
std::unique_ptr<capd::threading::ThreadPool> RosslerPMPool::threadPool;

//########################################################################
/**
 * A general algorithm that checks if P(X) satisfies the property 'condition'.
 * We assume that the domain X is covered by a grid of elements (newTasks) and 
 * to each element 'u' of this grid a task is assigned, which computes P(u) and/or derivatives of P.
 * Then we check is condition(P(u)) is satisfied. 
 * On failure, the domain 'u' is subdivided until subdivision limit (loop) is exceeded.  
 */ 
template<class Task, class Condition>
bool checkCondition(std::vector<Task*>& newTasks, Condition condition, int& counter)
{
  std::vector<Task*> completedTasks;

  bool result = false;
  int loop = 20; // max depth of subdivisions

  std::cout << "Initial number of tasks: " << newTasks.size() << std::endl << std::flush;
  while(loop and !result){
    for(Task* task : newTasks) RosslerPMPool::threadPool->process(task); // run all tasks
    for(Task* task : newTasks) task->join();             // join all tasks

    result = true;
    std::vector<Task*> tasks;

    for(Task* task : newTasks ) {
      bool r = condition(task); // check condition for completed task
      result = result and r;
      if(r)                     // On success, move it to completed tasks
        completedTasks.push_back(task);
      else{                     // otherwise split it (subdivide in (x,a) coordiantes)
        task->splitTask(tasks);
        delete task;
      }
    }
    std::swap(newTasks,tasks);
    std::cout << "New Boxes in subdivision: " << newTasks.size() << std::endl << std::flush;
    loop--;
  }

  for(Task* task : newTasks) delete task;

  std::swap(completedTasks,newTasks);
  std::cout << "Final number of boxes in subdivision: " << newTasks.size() << std::endl << std::flush;
  counter += newTasks.size();
  return result;
}
#endif
