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
   virtual double ComputePolyError();
   virtual double L1_error_projected_sol();
   virtual double L2_error_projected_sol();
   virtual double L1_error_average_sol();
   virtual double L2_error_average_sol();
   virtual double L1_error_nodal_gausspts();
   virtual double L2_error_nodal_gausspts();
   virtual double L1_error_nodal_gausspts_proj();
   virtual double L2_error_nodal_gausspts_proj();
   virtual double Compute_waveEnergy(double **in_Qn_);

   virtual void print_cont_vertex_sol();
   virtual void print_average_sol();
   virtual void dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                    ,double& L1_aver_sol_,double& L2_aver_sol_
                    ,double& L1_nodal_gausspts, double& L2_nodal_gausspts);

   virtual void dump_discont_sol();
   virtual void dump_timeaccurate_sol();
   virtual void dump_timeaccurate_errors();
   virtual void dump_timeaccurate_waveenergy(const double& in_E_ex_,
                                             const double& in_E_,
                                             const double& in_GG_ex_,
                                             const double& in_GG_);

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

   virtual double TimeAccurateExactSol_legendre_proj(const int &eID,
                                                const int &basis_k,
                                                const GaussQuad &quad_){}
   virtual void Compute_TimeAccurate_exact_sol(){}

private:
   double eval_basis_poly_derivative(const double& xi_pt, const int& basis_k_);
   double eval_local_du(const int eID, const double* q_, const double& xi_pt);
   double eval_local_du_fast(const int eID, const double* q_
                             , const int& position_);
   double Compute_common_du_flux(const double& dul, const double& dur);
   double eval_local_du_fluxproj(const int eID, const double* q_
                                 , const int& basis_k_);
   double eval_local_du_fluxproj_exact(const int eID, const double* q_
                                 , const int& basis_k_);
   double Compute_common_sol_jump(const double ul_, const double ur_);

   // CGR related functions:
   void load_cgr_scheme_matrices();
   void UpdateResidOneCell_cgr(const int &eID_,const double *qm1_,const double *q0_
                               ,const double *qp1_,double *resid_);

   double L1_norm_modal(double **Qn_, const int in_Nelem_
                        , const int in_ndof_);
   double L2_norm_modal(double **Qn_, const int in_Nelem_
                        , const int in_ndof_);

   double L1_error_modal(){}
   double L2_error_modal();

   void dump_gaussian_coeffs();

   void Reset_solver();

protected:

   double *u_sol_jump=nullptr;
   double C_lift=1.0;    // C lifting term multiplier for LDG, BR1, BR2/SIP
   double C1_lift=1.0;   // C1 lifting term multiplier for BR1
   double L1norm_modal_init_=0.0;
   double L2norm_modal_init_=0.0;

   double eta_penalty=1.0;  // penalty parameter
   GaussQuad quad_viscF_;
   int Nquad_viscFlux_=1;

   double **dLk=nullptr; //Derivatives of legendre polynomials at(-1), (1)

   typedef std::vector<double> vector;
   typedef std::vector<std::vector<double>> matrix;
   typedef std::vector<int> ivector;
   typedef std::vector<std::vector<int>> imatrix;

   matrix K0_s,Km1_s,Kp1_s;


};

#endif
