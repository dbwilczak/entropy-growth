#ifndef __CONJUGACY_DATA__
#define __CONJUGACY_DATA__

#include <vector>  
#include "TrappingRegion.h"
/**
 * This is an auxiliary class that collects non-rigorous algorithms for finding approximate extrema of the function 
 * f_z(x) = \pi_xP(x,0,z)
 * These extrema are then used to guess transition matrix of symbolic dynamics.
 * 
 * All algorithms take additional argument 'id' - which is the id of a thread from the thread pool.
 */   
struct Extrema{

  /// nonorigorous, compute x-th component of the Poincare map, that is 
  /// f_z(x) = \pi_xP(x,0,z)
  /// @returns f_z(x)
  static double f(double x, double z, int id){
    return RosslerPMPool::dRosslerPMPool->getSolver(id).pm({ x, 0., z})[0];
  }

  /// nonrigorous, compute f_z(x) and f_z'(x)
  /// @param [out] der = f_z'(x)
  /// @returns f_z(x)  
  static double df(double x, double z, double& der, int id){
    auto& pm = RosslerPMPool::dRosslerPMPool->getSolver(id).pm;
    DVector u({ x, 0., z});
    DMatrix D(3,3);
    u = pm(u,D);
    der = pm.computeDP(u,D)(1,1);
    return u[0];
  }
  
  /// Use simple bisection algorithm to find extremum of f_z in [x_1,x_2]
  /// (find zero of the derivative).
  /// It is assumed that f_z'(x) changes sign in the interval [x1,x2]
  /// @param [in] [x1,x2] - an interval that contains zero of the derivative of f_z
  /// @param [in] D1 = f_z'(x1), D2 = f_z'(x2) 
  static double bisect(double x1, double x2, double z, double D1, double D2, int id){
    double D, x;
    for(int i=0;i<20;++i) {
      x = (x1+x2)/2;
      df(x,z,D,id);
      if(D*D1<0){
        x2 = x;
        D2 = D;
      } else if(D*D2<0){
        x1 = x;
        D1 = D;
      } else
        return x;
    }
    return x;
  }

  /// The main algorithm: localize all extrema of f in the trapping region
  /// @param [in] a - parameter of the Rossler system
  /// @param [out] e - extrema found
  static void findExtrema(double a, double z, std::vector<double>& e, int id){
    auto& pm = RosslerPMPool::dRosslerPMPool->getSolver(id).pm;
    pm.getVectorField().setParameter("a",a);
    int N = 400;
    double D, prevD = 1;
    double prevX = TrappingRegion::xMax;
    double delta = (TrappingRegion::xMax-TrappingRegion::xMin)/N;
    for(int i=N; i>=0; i--){
      double x = TrappingRegion::xMin + i*delta;
      df(x,z,D,id);
      if(prevD*D<0) 
        e.push_back( bisect(x,prevX,z,D,prevD,id) );
      prevX = x;
      prevD = D;
    }
  }

  /// Auxiliary algoritm. It finds approximate bound of relevant extrema
  /// @param [in] a - range of parameters of the Rossler system
  /// @returns approximate location of extrema of f_Z(x) for all parameters from the interval a.
  static std::vector<interval> estimateRelevantExtrema(capd::interval a, int id){
    std::vector<double> tmp[6];
    double pts[][2] = {
        {a.rightBound(),TrappingRegion::zMin},{a.rightBound(),TrappingRegion::zMax},{a.rightBound(),TrappingRegion::zMid},
        {a.leftBound(),TrappingRegion::zMin},{a.leftBound(),TrappingRegion::zMax},{a.leftBound(),TrappingRegion::zMid}
      };
    for(int i=0;i<6;++i){
      findExtrema(pts[i][0],pts[i][1],tmp[i],id);
      std::sort(tmp[i].rbegin(),tmp[i].rend());
      double fx = f(tmp[i][0],pts[i][1],id);
      while(tmp[i].back()<=fx) tmp[i].pop_back();
    }
    // find minimal number of extrema
    int noProperExtrema = tmp[0].size();
    for(int i=1;i<6;++i){
      if(tmp[i].size()<noProperExtrema){
        for(int j=0;j<i;++j)
          tmp[j].pop_back();
        noProperExtrema--;
      }
    }
    // take hull
    std::vector<interval> result(tmp[0].begin(),tmp[0].end());
    for(int i=1;i<6;++i)
      for(unsigned j=0;j<result.size();++j)
        result[j] = intervalHull(result[j],tmp[i][j]);
  
    return result;
  }
};

/**
 * This class provides algorithms for finding possible conjugacy matrix for symbolic dynamics.
 * First, we find approximate minima and maxima of the function 
 * f(x) = \pi_xP(x,0,z)
 * for various parameter values 'a' of the Rossler system and for z=zMin,zMid,zMax.
 * 
 * Then, based on this data, it constructs h-sets for covering relations and the transition matrix for the shift dynamics. 
 * 
 * Such conjugacy is then validated in separate task.
 * 
 * The only public method and the main entry point is 'constructConjugacy'.
 */   
class ConjugacyData {
public:
  /// h-sets to be constructed
  std::vector< std::pair<double,double> > symbols;
  /// extrema of f_z(x)
  std::vector<double> extrema;
  /// computed transition matrix
  ZMatrix conjugacy;
  
  double minX, maxX; // subinterval of TrappingRegion::x-coordinate on which symbolic dynamics is constructed. 
  const capd::interval a;
  const double margin;
  int id; // thread id

  ConjugacyData(capd::interval a, double m) : a(a), margin(m) {}

  /// Main entry point. 
  void constructConjugacy(){
    // construct h-sets for symbolic dynamics
    this->findHSets();
    // based on this data, predict conjugacy matrix
    this->computeConjugacyMatrix();
    //~ for(auto set:  symbols) LOGGER(set);
  }

private:
  /// This is the most important routine. 
  /// It implements several heurestics for finding optimal location of h-sets for symbolic dynamics.
  /// 1. Localize estrema. This is easy part. 
  /// 2. The most difficult part is to adjust proper location of the most left and right sets. 
  void findHSets(){
    symbols.clear();
    extrema.clear();
    
    // 1. Find extrema
    auto extr = Extrema::estimateRelevantExtrema(a,id);
    for(auto e : extr) 
      extrema.push_back(e.rightBound());

    // 2. Find min-max of the attractor
    findMinMaxOfAttractor(extr);
    if(extr.back()<minX){
      extr.pop_back();
      extrema.pop_back();
    }

    // 3. Construct the most right set
    constructRightSet(extr);
        
    // 4. Construct intermediate sets from local extrema
    for(unsigned i=0;i<extrema.size()-1;++i)
      symbols.push_back(std::pair(extrema[i+1],extrema[i]));

    // 5. Construct the most left
    constructLeftSet(extr);
  }
  
  std::pair<double,double> getMinMax(double x){
    double fMin = TrappingRegion::xMax;
    double fMax = TrappingRegion::xMin;
    for(auto _a : {a.leftBound(),a.rightBound()}){
      for(auto z : {TrappingRegion::zMax,TrappingRegion::zMid,TrappingRegion::zMin}){ 
        RosslerPMPool::dRosslerPMPool->getSolver(id).vectorField.setParameter("a",_a);
        double fx = Extrema::f(x,z,id);
        fMin = capd::min(fMin,fx);
        fMax = capd::max(fMax,fx);
      }
    }
    return std::pair(fMin,fMax);
  }

  void constructLeftSet(std::vector<interval>& extr){
    unsigned noExtrema = extrema.size();
    double x = minX + 0.1*margin*(extr.back().leftBound()-minX);
    double fx = getMinMax(x).first;
    if( noExtrema%2==1 )
    {
      // check if the set covers at least itself
      if(fx>symbols.back().first + margin*(symbols.back().first-x)){
        //~ std::cout << "case 1:\n";
        symbols.push_back(std::pair(x,extrema.back()));
      } else { // last set is just after minimum and seems to be stable
        //~ std::cout << "case 2:\n";
        extrema.pop_back();
      }
    } else {
      // we are after maximum. Check if the set covers the right set
      if(fx<symbols[0].first - margin*(symbols[0].second-symbols[0].first)){
        //~ std::cout << "case 3:\n";
        symbols.push_back(std::pair(x,extrema.back()));        
      } else {
        //~ std::cout << "case 4:\n";
          extrema.pop_back();
      }
    }
  }

  // There are two cases for the most right set
  void constructRightSet(std::vector<interval>& extr){
    // Case 1. 
    // - number of extrema is even and the first set will always be covered by the last.
    // - the most left set has no chance to cover the right set
    // then the size of right set is irrelevant.
    double x = minX + margin*(extr.back().leftBound()-minX);
    double fL = getMinMax(x).first;
    if(extr.size()%2==0 or fL<extr[0]-margin*(extrema[0]-extrema[1]))
    {
      double R = maxX - margin*(maxX-extrema[0]);
      symbols.push_back(std::pair(extrema[0],R));
      return;
    }
    
    // 2. Now the number of extrema is odd and there is possibility for the left set to cover right set.
    // In this case we want to make the set as small as possible but to keep the chance of covering as much as possible sets.
    double R = maxX;
    double fR = getMinMax(R).first;
    
    unsigned i = 0;
    while( !(fR>extr[i]+0.01) )  i++;
    double L = extrema[0];
    for(int j=0;j<20;++j){
      double x = (L+R)/2;
      double fx = getMinMax(x).first;
      if( !(fx>extr[i]+0.01) )
        L = x;
      else
        R = x;
    }
    symbols.push_back(std::pair(extrema[0],R));
  }

  void findMinMaxOfAttractor(std::vector<interval>& extr){
    unsigned noExtrema = extrema.size();
    minX = getMinMax(extr[0].leftBound()).second;
    minX = capd::max(minX,getMinMax(extr[0].rightBound()).second);
    for(unsigned i=2;i<noExtrema;i+=2) {
      minX = capd::max(minX,getMinMax(extr[i].leftBound()).second);
      minX = capd::max(minX,getMinMax(extr[i].rightBound()).second);
    }
    // There are two cases for maxX. 
    // 2a) number of extrema>1. Then the local maxima give the bound.
    // 2b) there is only one minimum. Then we have to take into account value at xMin and xMax to obtain maximum.
    if(noExtrema>1){
      maxX = getMinMax(extr[1].leftBound()).first;
      maxX = capd::min(maxX,getMinMax(extr[1].rightBound()).first);
      for(unsigned i=3;i<noExtrema;i+=2) {
        maxX = capd::min(maxX,getMinMax(extr[i].leftBound()).first);
        maxX = capd::min(maxX,getMinMax(extr[i].rightBound()).first);
      }
    } else {
      maxX = getMinMax(minX).first;
      maxX = capd::max(getMinMax(TrappingRegion::xMax).first,maxX);
    }
  }
  
  void computeConjugacyMatrix(){
    unsigned N = symbols.size();
    conjugacy = ZMatrix(N,N);
    unsigned i;
    for(i=0;i<N;++i){
      auto x1 = getMinMax(symbols[i].first);
      auto x2 = getMinMax(symbols[i].second);
      for(unsigned j=0;j<N;++j){
        if(x1.second<symbols[j].first and x2.first>symbols[j].second) 
          conjugacy(i+1,j+1) = 1;
        if(x2.second<symbols[j].first and x1.first>symbols[j].second) 
          conjugacy(i+1,j+1) = -1;
      }
    }
    bool emptyColumn = true;
    for(i=1;i<=N and emptyColumn;i++)
      if(conjugacy(i,N)!=0) emptyColumn = false;
    if(emptyColumn){
      symbols.pop_back();
      ZMatrix c(N-1,N-1);
      for(i=1;i<N;++i)
        for(unsigned j=1;j<N;++j)
          c(i,j) = conjugacy(i,j);
          conjugacy = c;
    }
  }
  
};

#endif
