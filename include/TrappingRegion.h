#ifndef ROSSLER_TEMPLATES_TRAPPING_REGION_H
#define ROSSLER_TEMPLATES_TRAPPING_REGION_H

#include "capd/capdlib.h"
using namespace capd;

/**
 * This is a data structure, which for a given parameter (possible interval) 'a' represents 
 * a trapping region for the Poincare map. 
 * The trapping region is a rectangle alligned to axes.
 * Right bottom and right left corners are fixed for all parameter values --- see constants above.
 * The 'x'-coordinate of left edge is computed as:
 * g_Lmax - (g_Lmax-g_Lmin)*(a-g_amin)/(g_amax-g_amin) 
 */ 
struct TrappingRegion{
  // range of parameter we examine
  constexpr static double aMin = 0.12;
  constexpr static double aMax = 0.3659;

  // coordinates of a trapping region
  constexpr static double zMin = 0.004;
  constexpr static double zMax = 0.011;
  constexpr static double zMid = (zMin+zMax)/2;
  constexpr static double xMin = -30.53;
  constexpr static double xMax = -3;  

  std::deque<double> grid;
  capd::interval a;
  
  explicit TrappingRegion(capd::interval a){
    this->a = a;
    double x=xMax;
    grid.push_back(x);
    while(x>=xMin){
      x -=  exp(x/5);
      grid.push_front(x>xMin ? x : xMin);
    }
  }
  
  void setParameter(interval a){
    this->a = a;
  }
  
  capd::IVector gridElement(int k) const {
    return capd::IVector({ capd::interval(grid[k],grid[k+1]), capd::interval(0.), capd::interval(zMin,zMax), a });
  }

  capd::IVector getTrappingRegion() const {
    return capd::IVector({ capd::interval(xMin,xMax), capd::interval(0.), capd::interval(zMin,zMax), a });    
  }

  capd::IVector getLeftEdge() const {
    return capd::IVector( { capd::interval(xMin), capd::interval(0.), capd::interval(zMin,zMax), a } );    
  }
  
  capd::IVector getRightEdge() const {
    return capd::IVector( { capd::interval(xMax), capd::interval(0.), capd::interval(zMin,zMax), a } );    
  }

  capd::IVector getTopEdge() const {
    return capd::IVector({ capd::interval(xMin,xMax), capd::interval(0.), capd::interval(zMax), a });    
  }

  capd::IVector getBottomEdge() const {
    return capd::IVector({ capd::interval(xMin,xMax), capd::interval(0.), capd::interval(zMin), a });    
  }
};

#endif
