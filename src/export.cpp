#include "main.h"
#include "ConjugacyData.h"

using namespace std;
using namespace capd;

void gridEdges(std::vector<IVector>& grid, double a){
  TrappingRegion tr(a);
  for(unsigned i=1;i<tr.grid.size();i++){
    IVector u = tr.gridElement(i-1); // a set of the form [x]\times [0] \times [z]
    IVector u1=u, u2=u;
    u1[2] = u1[2].rightBound(); // u1 = [x]\times [0] \times [zMaxz]
    u2[2] = u2[2].leftBound();  // u1 = [x]\times [0] \times [zMin]
    grid.push_back(u1);
    grid.push_back(u2);
  }
  grid.push_back(tr.getLeftEdge());
  grid.push_back(tr.getRightEdge());
}

//########################################################################
void exportEnclosureData(double a, const char* filename){
  RosslerPMPool::iRosslerPMPool->setOrder(20);
  RosslerPMPool::iRosslerPMPool->setTolerance(1e-21,1e-21);
  std::vector<IVector> grid, tmp;
  gridEdges(grid,a);

  auto cond = [](const IVector& u) { 
      return 
             u[0].leftBound()  >= TrappingRegion::xMin and 
             u[0].rightBound() <= TrappingRegion::xMax and
             u[2].leftBound()  >= TrappingRegion::zMin and 
             u[2].rightBound() <= TrappingRegion::zMax;    
  };

  FILE* out = fopen(filename,"w");
  while(grid.size()>0){
    for(int i=0;i<grid.size();++i){
      C0HOTripletonSet s(grid[i]);
      IVector v = RosslerPMPool::iRosslerPMPool->getSolver(0).pm(s);
      if( cond(v) ){
        fprintf(out,"%f %f %f %f\n",v[0].leftBound(),v[2].leftBound(),v[0].rightBound(),v[2].rightBound());
      } else {
        IVector u1 = grid[i];
        IVector u2 = grid[i];
        double x = (grid[i][0].leftBound()+grid[i][0].rightBound())/2;
        u1[0].setRightBound(x);
        u1[0].setLeftBound(x);
        tmp.push_back(u1);
        tmp.push_back(u2);
      }
    }
    grid.clear();
    std::swap(tmp,grid);
    LOGGER(grid.size());
  }
  
  fclose(out);
}

//########################################################################
void exportEnclosureData(){
  exportEnclosureData(TrappingRegion::aMin,"./data/aminEnclosure.dat");
  //~ exportEnclosureData(0.3,"./data/Enclosure03.dat");
  exportEnclosureData(0.325,"./data/Enclosure0325.dat");
}


void evalEdge(double x, double a, FILE* out, int N){
  static int i = 0;
  interval d = (interval(TrappingRegion::zMax)-interval(TrappingRegion::zMin))/N;
  IVector v;
  for(int k=0;k<N;++k){
    IVector u{x+interval(-1,1)*pow(2.,-52),0.,TrappingRegion::zMin + interval(k,k+1)*d,a};
    C0HOTripleton s(u);
    if(k==0)
      v = RosslerPMPool::iRosslerPMPool->getSolver(0).pm(s);
    else 
      v = intervalHull(v,RosslerPMPool::iRosslerPMPool->getSolver(0).pm(s));
  }
  fprintf(out,"%.17f %.17f %.17f %.17f\n",v[0].leftBound(),v[2].leftBound(),v[0].rightBound(),v[2].rightBound());
  cout <<" X_{" << i << "} & \\subset & " << v[0] << ", \\\\\n";
  ++i;
}

//########################################################################
void exportCovrel(double a, const char* filename, int N = 40){
  ConjugacyData c(a,1e-3);
  RosslerPMPool::iRosslerPMPool->setOrder(20);
  RosslerPMPool::iRosslerPMPool->setTolerance(1e-21,1e-21);
  c.id = 0;
  c.constructConjugacy();
  FILE* out = fopen(filename,"w");
  for(int i=0;i<c.symbols.size();++i){
    evalEdge(c.symbols[i].second,a,out,N);
    //~ cout << c.symbols[i].first << " " << c.symbols[i].second << endl;
  }
  
  evalEdge(c.symbols.back().first,a,out,N);
  fclose(out);
  for(int i=c.symbols.size()-1;i>=0;--i){
    cout << c.symbols[i].first << ",";
  }
  cout << c.symbols[0].second;
  
}

//########################################################################

void exportCurve(double a, const char* filename){
  RosslerPMPool::dRosslerPMPool->setParameter("a",a);
  int N = 7000;
  
  FILE* out = fopen(filename,"w");
  for(int i=0;i<=N;++i){
    double x = TrappingRegion::xMin + i*(TrappingRegion::xMax-TrappingRegion::xMin)/N;
    fprintf(out,"%f %f\n",x,Extrema::f(x,TrappingRegion::zMid,0));
  }
  fclose(out);
}

void eigenvalues(double a, double b, double c, double d){
  double s = sqrt(a*a + 4*b*c - 2*a*d + d*d);
  double lambda1 = 0.5* (a + d - s);
  double lambda2 = 0.5* (a + d + s);
  LOGGER(lambda1);
  LOGGER(lambda2);
}

void findBifPoint(){
  auto& PM = *RosslerPMPool::DPMapFactory("par:b,c;var:x,y,z,a;fun:-(y+z),x+a*y,b+z*(x-c),0;",20,2,capd::DCoordinateSection(4,1),capd::poincare::PlusMinus).createSolver();
  auto& pm = PM.pm;
  PM.vectorField.setParameter("b",0.2);
  PM.vectorField.setParameter("c",15);
  PM.solver.setOrder(10);
  
  //~ DVector u{-24.98615641824929,0,0.005003695953767296,0.2445890212249042};
  //~ DVector u{-27.90101006311546,0,0.004663547021688063,0.311986650909318};
  //~ DVector u{-29.22211573599117,0,0.004524155280447693,0.3405236989996505};
  //~ DVector u  {-29.89377365336397,0,0.004456435009829962,0.3546641965836549};
  DVector u {-30.25817997684925,0,0.004420535075303823,0.3623180183723328};
  DMatrix Dflow(4,4), D(4,4);
  DHessian Hflow(4,4), H(4,4);
  for(int i=0;i<50;++i){
    DVector v = pm(u,Dflow,Hflow);
    pm.computeDP(v,Dflow,Hflow,D,H);
    double a = 2*H(0,0,0)*(D(3,3)-1) + (D(1,1)-1)*H(2,0,2) - H(0,0,2)*D(3,1) - 2*D(1,3)*H(2,0,0);
    double b = 2*H(0,2,2)*(D(1,1)-1) + (D(3,3)-1)*H(2,0,2) - H(0,0,2)*D(1,3) - 2*D(3,1)*H(2,3,3);
    double c = H(0,1,3)*(D(3,3)-1) + (D(1,1)-1)*H(2,2,3) - H(0,2,3)*D(3,1) - D(1,3)*H(2,0,3);
    DMatrix A{
        {D(1,1)-1,D(1,3),D(1,4)},
        {D(3,1),D(3,3)-1,D(3,4)},
        {a,b,c}
      };
    DVector defect{v[0]-u[0],v[2]-u[2], (D(1,1)-1)*(D(3,3)-1)- D(1,3)*D(3,1)};
    LOGGER(defect);
    DVector e = matrixAlgorithms::gauss(A,defect);
    LOGGER(e);
    u[0] -= e[0];
    u[2] -= e[1];
    u[3] -= e[2];
  }
  eigenvalues(D(1,1),D(1,3),D(3,1),D(3,3));
  LOGGER(u);
}

/*
    std::ofstream out;
    out.open("./data/lbranch1.dat");
    out.precision(17);
    for(auto* t : lBranch)
      out << t->u[0].mid().leftBound() << " " << t->u[3].mid().leftBound()<< endl;
    out.close();
    out.open("./data/rbranch1.dat");
    out.precision(17);
    for(auto* t : rBranch)
      out << t->u[0].mid().leftBound() << " " << t->u[3].mid().leftBound()<< endl;
    out.close();  

*/

double getEntropy(ZMatrix _A){
  int N = _A.numberOfRows();
  DMatrix A(_A);
  DVector u(N);
  for(int i=0;i<N;++i) u[i] = i;
  double L = 0;
  for(int i=0;i<10000;++i){
    u.normalize();
    u = A*u;
    L = capd::max(L,u.euclNorm());
  }
  return log(L);
}
  

void readMatrices(){
  ifstream in("./out/macierze.dat");
// ConjugacyResult is a list of pairs (interval,transition matrix)
  typedef std::vector< std::pair<interval,ZMatrix> > ConjugacyResult;
  ConjugacyResult cr;
  interval a;
  ZMatrix M;
  while(!in.eof()){
    in >> a;
    if(isSingular(a)) break;
    in >> M;
    cr.push_back(std::pair(a,M));
  }
  in.close();
  cout << cr.size() << endl;
  cout.precision(6);
  const int N = 3; // print 3 data in one row
  const int n = 1+cr.size()/N;
  cout << n << endl;
  for(int i=0;i<n;++i){
    for(int j=0;j<N;++j){
      int k = i+n*j;
      if(k<cr.size()){
        auto [a,M] = cr[k];
        std::cout << "$" << a.leftBound() << "$ & $" << M.numberOfRows() << "$ & $" << getEntropy(M) << "$";
      }
      if(j==N-1){
        std::cout << "\\\\ \\hline\n";
      } else {
        std::cout << " & ";
      }
    }
  }
}

void printEntropy(){
  ifstream in("./data/ev.dat");
  DMatrix M;
  in >> M;
  in.close();
  cout << M;
  
  cout.precision(6);
  const int N = 3; // print 3 data in one row
  const int n = 1+M.numberOfRows()/N;
  cout << n << endl;
  for(int i=0;i<n;++i){
    for(int j=0;j<N;++j){
      int k = i+n*j+1;
      if(k<=M.numberOfRows()){
        std::cout << "$" << k << "$ & $" << M(k,1) << "$ & $" << M(k,2) << "$ & $" << M(k,3) << "$";
      }
      if(j==N-1){
        std::cout << "\\\\ \\hline\n";
      } else {
        std::cout << " & ";
      }
    }
  }

}

int main(int argc, char* argv[])
{
  cout.precision(16); 
  try{
    // Init a pool of objects that compute Poincare map for the Rossler system.
    // The size of the pool=hardware_concurrency is the maximal available number of threads on the machine.
    RosslerPMPool::init(std::thread::hardware_concurrency());
    printEntropy();
   
    //~ findFixedPoint();
    //~ findBifPoint();


    //~ exportCurve(0.11,"data/curve.dat");

    // exporting data for pictures
    //~ exportEnclosureData();
    //~ exportCovrel(0.3659,"data/aminCovrel.dat",40);

    //~ exportCurve(TrappingRegion::aMin,"data/aminCurve.dat");
    //~ exportCurve(0.2,"data/curve2.dat");
    //~ exportCurve(0.27,"data/curve27.dat");
    //~ exportCurve(0.3,"data/curve3.dat");
    //~ exportCurve(0.32,"data/curve32.dat");
    //~ exportCurve(0.333,"data/curve333.dat");
    //~ exportCurve(0.35,"data/curve35.dat");
    //~ exportCurve(0.356,"data/curve356.dat");
    //~ exportCurve(0.366,"data/curve366.dat");
    //~ exportCurve(0.3658,"data/curve3658.dat");
    //~ exportCurve(0.3659,"data/curve3659.dat");
    //~ exportCurve(0.3666,"data/curve3666.dat");
    //~ exportCurve(TrappingRegion::aMax,"data/amaxCurve.dat");

  }catch(const exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }

  return 0; 
}
