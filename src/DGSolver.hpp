#ifndef DGSOLVER_H
#define DGSOLVER_H


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

   void setup_solver(GridData& meshdata_, SimData& simdata_);

   void InitSol();
   void UpdateResid(double **Resid_, double **Qn_);
   void UpdateSolution(double **Qn_);
   void ComputeError();
   void Compute_vertex_sol();
   void Compute_exact_sol();
   void UpdatePhyTime(const double& dt_);

   void SetPhyTime(const double& time_);
   double GetPhyTime();
   double GetTimeStep();
   double GetCFL();
   unsigned int GetNdof();
   double** GetNumSolution();
   double* GetVertexNumSol();
   double* GetExactSolution();

   void print_num_vertex_sol();
   void print_exact_sol();
   void print_exact_average_sol();
   void print_num_average_sol();

protected:
   void UpdateResidOneCell(const int& cellid, double* q_, double* resid_);

   double Compute_common_flux(const double& ql, const double& qr,
                               const double& wave_speed
                               , const double& upwind_Beta_);
   void Compute_flux_upw();
   void get_left_right_sol();
   void Rusanov_flux();
   void Roe_flux();

   double eval_init_sol(const double& xx);
   double eval_basis_poly(const double& xi_, const int& basis_k_);
   double eval_basis_norm_squared(const int& basis_k_);
   double evalSolution(const double* q_, const double& xi_pt);

   double eval_localflux_proj(const double* q_, const unsigned int& basis_k_);
   double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_);

   void CalcTimeStep();
   void CalcLocalTimeStep();

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
   unsigned int Ndof = 1;

   double **Qn=nullptr;      // Nelem * Ndof long

   double *Q_exact=nullptr;  // Nfaces long

   double *Qv=nullptr;       // Nfaces long

   double *flux_com=nullptr;  // common interface flux, Nfaces long

   //double **Resid=nullptr;

   //double *Mnn=nullptr;

   //double *phy_n=nullptr;

   double phy_time=0.0;
   double time_step=1e-5;
   double CFL=1.0;

   //int Beta_=1;  // upwind parameter

};

#endif
