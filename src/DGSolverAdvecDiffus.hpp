#ifndef DGSOLVERADVECDIFFUS_H
#define DGSOLVERADVECDIFFUS_H


#include "DGSolver.hpp"

//struct SimData;
//struct GridData;

class DGSolverAdvecDiffus:public DGSolver{

public:

//  Construction Functions :
   DGSolverAdvecDiffus(void){}
   virtual ~DGSolverAdvecDiffus(void);

   virtual void setup_solver(GridData& meshdata_, SimData& simdata_);

   virtual void InitSol();
   virtual void UpdateResid(double **Resid_, double **Qn_);

   virtual void Compute_vertex_sol();
   virtual double ComputePolyError();
   virtual double L1_error_projected_sol();
   virtual double L2_error_projected_sol();
   virtual double L1_error_average_sol();
   virtual double L2_error_average_sol();
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
   virtual void dump_timeaccurate_sol();
   virtual void dump_timeaccurate_errors(){}

   virtual void UpdateResidOneCell(const int& cellid, double* q_
                           , double* resid_);
   virtual double eval_init_sol(const double& xx);
   virtual double eval_exact_sol(double& xx);
   virtual double eval_basis_poly(const double& xi_, const int& basis_k_);
   virtual double eval_basis_norm_squared(const int& basis_k_);
   virtual void setup_basis_interpolation_matrices();

   virtual double evalSolution(const double* q_, const double& xi_pt);
   virtual double evalSolution(const double *q_, const int& position_);
   virtual double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_);
   virtual double ExactSol_legendre_proj(const int &eID,
                                 const int &basis_k,
                                  const GaussQuad &quad_);
   virtual void CalcTimeStep();
   virtual void ComputeExactSolShift();
   virtual void Compute_exact_vertex_sol();
   virtual void Compute_projected_exact_sol();

   virtual void compute_uniform_cont_sol();
   virtual double compute_totalVariation();

protected:

   double eval_basis_poly_derivative(const double& xi_pt, const int& basis_k_);
   double eval_local_du(const int eID, const double* q_, const double& xi_pt);
   double eval_local_du_fast(const int eID, const double* q_
                             , const int& position_);
   double Compute_common_du_flux(const double& dul, const double& dur);
   double eval_local_du_fluxproj(const int eID, const double* q_
                                 , const int& basis_k_);
   double eval_local_du_fluxproj_exact(const int eID, const double* q_
                                 , const int& basis_k_);
   double Compute_common_sol_jump(const double& ul_, const double& ur_);

   // Inviscid Flux functions:
   double Compute_common_invflux(const double& ql, const double& qr,
                               const double& wave_speed
                               , const double& upwind_Beta_);
   double eval_local_invflux_proj(const double* q_
                                 , const int& basis_k_);
   double eval_local_invflux_proj_exact(const double* q_
                                 , const int& basis_k_);
   double Rusanov(const double &Ql, const double &Qr);
   double eval_burgers_invflux(const double& xi_pt, const double *q_);
   double eval_init_u_decay_burger_turb(const double& xx_);

   void Reset_solver();

protected:

   double *u_sol_jump=nullptr;
   double r_lift=1.0;
   double e_penalty=1.0;
   double *viscflux_com=nullptr;
   GaussQuad quad_invF_;
   GaussQuad quad_viscF_;
   int Nquad_invFlux_=1;
   int Nquad_viscFlux_=1;

   //double *Q_cont_sol=nullptr;

   double **dLk=nullptr; //Derivatives of legendre polynomials at(-1), (1)
};

#endif
