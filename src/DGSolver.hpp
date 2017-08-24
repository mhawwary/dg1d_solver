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

    double T_period=1;

//  Construction Functions :
   DGSolver(){}
   virtual ~DGSolver(){}

   virtual void setup_solver(GridData& meshdata_, SimData& simdata_)=0;

   virtual void InitSol()=0;
   virtual void UpdateResid(double **Resid_, double **Qn_)=0;
//   virtual void UpdateSolution(double **Qn_)=0;

   virtual void Compute_vertex_sol()=0;
   virtual void Compute_cont_sol()=0;
   virtual double ComputePolyError()=0;
   virtual double L1_error_projected_sol()=0;
   virtual double L2_error_projected_sol()=0;
   virtual double L1_error_average_sol()=0;
   virtual double L2_error_average_sol()=0;
   virtual double L2_error_nodal_disc_sol()=0;
   virtual double L1_error_nodal_gausspts()=0;
   virtual double L2_error_nodal_gausspts()=0;
   virtual double L1_error_nodal_gausspts_proj()=0;
   virtual double L2_error_nodal_gausspts_proj()=0;

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

   virtual void print_cont_vertex_sol()=0;
   virtual void print_average_sol()=0;
   virtual void dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                    ,double& L1_aver_sol_,double& L2_aver_sol_
                    ,double& L1_nodal_gausspts, double& L2_nodal_gausspts)=0;
   virtual void dump_discont_sol()=0;
   virtual void dump_timeaccurate_sol()=0;

protected:

   virtual void UpdateResidOneCell(const int& cellid, double* q_
                           , double* resid_)=0;

//   virtual double Compute_common_flux(const double& ql, const double& qr,
//                               const double& wave_speed
//                               , const double& upwind_Beta_)=0;
//   virtual void Compute_flux_upw()=0;
//   virtual void get_left_right_sol()=0;
//   virtual void Rusanov_flux()=0;
//   virtual void Roe_flux()=0;

   virtual double eval_init_sol(const double& xx)=0;
   virtual double eval_exact_sol(double& xx)=0;
   virtual double eval_basis_poly(const double& xi_, const int& basis_k_)=0;
   virtual double eval_basis_norm_squared(const int& basis_k_)=0;
   virtual double evalSolution(const double* q_, const double& xi_pt)=0;

//   virtual double eval_localflux_proj(const double* q_
//                              , const int& basis_k_)=0;

   virtual double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_)=0;

   virtual double ExactSol_legendre_proj(const int &eID,
                                 const int &basis_k,
                                  const GaussQuad &quad_)=0;

   virtual void CalcTimeStep()=0;
//   virtual void CalcLocalTimeStep()=0;

   virtual void ComputeExactSolShift()=0;

   virtual void Compute_exact_vertex_sol()=0;
   virtual void Compute_projected_exact_sol()=0;

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

   int Nquad_=5; // Gauss Quadrature rules

};

#endif
