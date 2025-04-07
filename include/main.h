#ifndef ROSSLER_TEMPLATES_MAIN_H
#define ROSSLER_TEMPLATES_MAIN_H

#include <iostream>
#include <chrono>
#include "capd/dynset/C0HOSet.hpp"
#include "capd/dynset/C0TripletonSet.hpp"
#include "TrappingRegion.h"
#include "RosslerPMPool.h"

#define LOGGER(x) std::cout << std::boolalpha << (#x) << ": " << (x) << std::endl;

typedef capd::vectalg::Matrix<int,0,0> ZMatrix;
typedef capd::dynset::FactorReorganization< capd::dynset::PartialQRWithPivoting<3> > Policies;
typedef capd::dynset::C0TripletonSet<capd::IMatrix,Policies> C0Tripleton;
typedef capd::dynset::C0HOSet<C0Tripleton> C0HOTripleton;

inline std::chrono::time_point<std::chrono::system_clock> tic(const char message[]){
  auto start = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(start);
  std::cout << "\nStart " << message << ": " << std::ctime(&t) << std::endl;
  std::cout << "Number of threads: " << std::thread::hardware_concurrency() << "\n####################################\n\n";
  return start;
}

inline void tac(std::chrono::time_point<std::chrono::system_clock>& start){
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;  
  auto t = std::chrono::system_clock::to_time_t(end);
  std::cout << "Computation completed: " << std::ctime(&t)
            << "Elapsed time: " << elapsed_seconds.count()/60. << " minutes\n";    
}

#endif
