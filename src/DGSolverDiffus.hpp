#ifndef DGSOLVERDIFFUS_H
#define DGSOLVERDIFFUS_H


#include "DGSolver.hpp"

//struct SimData;
//struct GridData;

class DGSolverDiffus:public DGSolver{

public:

//  Construction Functions :
   DGSolverDiffus(void){}
   virtual ~DGSolverDiffus(void);

   virtual void setup_solver(GridData& meshdata_, SimData& simdata_);

   virtual void InitSol();
   virtual void UpdateResid(double **Resid_, double **Qn_);
//   virtual void UpdateSolution(double **Qn_);

   virtual void Compute_vertex_sol();
   virtual void Compute_cont_sol(){}
   virtual double ComputePolyError();
   virtual double L1_error_projected_sol();
   virtual double L2_error_projected_sol();
   virtual double L1_error_average_sol();
   virtual double L2_error_average_sol();
   virtual double L2_error_nodal_disc_sol();
   virtual double L1_error_nodal_gausspts();
   virtual double L2_error_nodal_gausspts();
   virtual double L1_error_nodal_gausspts_proj();
   virtual double L2_error_nodal_gausspts_proj();

   virtual void print_cont_vertex_sol();
   virtual void print_average_sol();
   virtual void dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                    ,double& L1_aver_sol_,double& L2_aver_sol_
                    ,double& L1_nodal_gausspts, double& L2_nodal_gausspts);
   virtual void dump_discont_sol();
   virtual void dump_timeaccurate_sol(){}

   virtual void UpdateResidOneCell(const int& cellid, double* q_
                           , double* resid_);

//   virtual void Compute_flux_upw();
//   virtual void get_left_right_sol();
//   virtual void Rusanov_flux();
//   virtual void Roe_flux();

   virtual double eval_init_sol(const double& xx);
   virtual double eval_exact_sol(double& xx);
   virtual double eval_basis_poly(const double& xi_, const int& basis_k_);
   virtual double eval_basis_norm_squared(const int& basis_k_);
   virtual double evalSolution(const double* q_, const double& xi_pt);

   virtual double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_);

   virtual double ExactSol_legendre_proj(const int &eID,
                                 const int &basis_k,
                                  const GaussQuad &quad_);

   virtual void CalcTimeStep();
//   virtual void CalcLocalTimeStep();

   virtual void ComputeExactSolShift();

   virtual void Compute_exact_vertex_sol();
   virtual void Compute_projected_exact_sol();

protected:

   double eval_basis_poly_derivative(const double& xi_pt, const int& basis_k_);

   double eval_local_du(const int eID, const double* q_, const double& xi_pt);

   double Compute_common_du_flux(const double& dul, const double& dur);

   double eval_local_du_fluxproj(const int eID, const double* q_
                                 , const int& basis_k_);

   double Compute_common_sol_jump(const double& ul_, const double& ur_);

   void Reset_solver();

protected:

   double *u_sol_jump=nullptr;
   double r_lift=1.0;
   double e_penalty=1.0;

};

#endif
