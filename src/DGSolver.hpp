﻿#ifndef DGSOLVER_H
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

    double T_period=1;

//  Construction Functions :
   DGSolver(void);
   ~DGSolver(void);

   void setup_solver(GridData& meshdata_, SimData& simdata_);

   void InitSol();
   void UpdateResid(double **Resid_, double **Qn_);
   void UpdateSolution(double **Qn_);

   void Compute_vertex_sol();
   double ComputePolyError();
   double Compute_projected_sol_error();
   double ComputeAverageError();
   double ComputeDiscNodalError();

   void UpdatePhyTime(const double& dt_){

       phy_time += dt_;

       return;
   }

   void SetPhyTime(const double &time_){

       phy_time=time_;

       return;
   }

   double GetTimeStep(){

       return time_step;
   }

   double GetLastTimeStep(){

       return last_time_step;
   }

   double GetCFL(){

       return CFL;
   }

   int GetNdof(){

       return Ndof;
   }

   double** GetNumSolution(){

       return Qn;
   }

   double* GetVertexNumSol(){

       return Qv;
   }

   double* GetExactSolution(){

       return Q_exact;
   }

   double GetPhyTime(){

       return phy_time;
   }

   void print_cont_vertex_sol();
   void print_average_sol();
   void dump_errors(double& proj_sol_L2, double &aver_L2);
   void dump_discont_sol();

protected:

   void UpdateResidOneCell(const int& cellid, double* q_
                           , double* resid_);

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

   double eval_localflux_proj(const double* q_
                              , const int& basis_k_);

   double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_);

   double ExactSol_legendre_proj(const int &eID,
                                 const int &basis_k,
                                  const GaussQuad &quad_);

   void CalcTimeStep();
   void CalcLocalTimeStep();

   void Reset_solver();
   void ComputeExactSolShift();

   void Compute_exact_vertex_sol();
   void Compute_projected_exact_sol();

protected:

   GridData *grid_=nullptr;
   SimData *simdata_=nullptr;

   double *xi=nullptr;

   int Ndof = 1;

   double **Qn=nullptr;      // Nelem * Ndof long

   double *Q_exact=nullptr;  // Nfaces long

   double **Qex_proj=nullptr; // projected exact solution , Nelem long

   double *Qv=nullptr;       // Nfaces long

   double *flux_com=nullptr;  // common interface flux, Nfaces long

   double phy_time=0.0;
   double time_step=1e-5;
   double last_time_step=1e-5;
   double CFL=1.0;

   double exact_sol_shift=0.;
   double wave_length_=0.;

};

#endif
