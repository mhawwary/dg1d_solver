#include"GirdData.h"
#include"SimData.hpp"
#include"quadrature.h"
#include"general_tools.h"
#include"global_var.h"

//struct SimData;
//struct GridData;

class DGSolver{

public:
//  Construction Functions :
   DGSolver(void);
   ~DGSolver(void);

   void setup_solver(GridData *meshdata_, SimData* simdata_);

   void InitSol();
   void UpdateResid(double **Resid_, double **Qn_);
   void UpdateSolution();
   void ComputeError();
   void Compute_vertex_sol();
   void Compute_exact_sol();

protected:
   void UpdateResidOneCell(const int& cellid, double* q_, double* resid_);
   void Compute_common_flux();
   void Compute_flux_upw();
   void get_left_right_sol();
   void Rusanov_flux();
   void Roe_flux();

   double eval_init_sol(const double& xx);
   double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_);
   void Reset_solver();


protected:

   GridData *grid_=nullptr;
   SimData *simdata_=nullptr;

   //int Nelem=1;
   //int Nfaces=1;

   //double *X=nullptr;
   //double *Xc=nullptr;
   //double *h_j=nullptr;

   double *xi=nullptr;

   //int poly_order=1;
   int Ndof = 1;

   double **Qn=nullptr;

   double *Q_exact=nullptr;

   double *flux_com=nullptr;  // common interface flux

   double **Resid=nullptr;

   double *Mnn=nullptr;

   double *phy_n=nullptr;

   //int Beta_=1;  // upwind parameter

};
