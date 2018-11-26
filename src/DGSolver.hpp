#ifndef DGSOLVER_H
#define DGSOLVER_H


#include"GridData.h"
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
   virtual void UpdateResid(double **Resid_, double **Qn_)=0; // updating residual
//   virtual void UpdateSolution(double **Qn_)=0;

   virtual void Compute_vertex_sol()=0; // compute num sol at element nodes and average
   virtual double ComputePolyError()=0;
   virtual double L1_error_projected_sol()=0; // L1 of projected exact and numerical solutions
   virtual double L2_error_projected_sol()=0; // L2 of projected exact and numerical solutions
   virtual double L1_error_average_sol()=0;   // L1 of the averages of both exact and numerical solutions
   virtual double L2_error_average_sol()=0;   // L2 of the averages of both exact and numerical solutions
   virtual double L1_error_nodal_gausspts()=0; // L1 at gauss points, but using non-projected exact solution
   virtual double L2_error_nodal_gausspts()=0; // L2 at gauss points, but using non-projected exact solution
   virtual double L1_error_nodal_gausspts_proj()=0; // L1 at gauss points, but using projected exact and numerical solutions
   virtual double L2_error_nodal_gausspts_proj()=0; // L2 at gauss points, but using projected exact and numerical solutions
   virtual double Compute_waveEnergy(double **in_Qn_)=0;

   inline void UpdatePhyTime(const double& dt_){  phy_time += dt_; }
   inline void SetPhyTime(const double &time_){ phy_time=time_; }
   inline double GetTimeStep(){  return time_step; }
   inline double GetLastTimeStep(){ return last_time_step; }
   inline double GetCFL(){  return CFL; }
   inline int GetNdof(){ return Ndof; }
   inline double** GetNumSolution(){ return Qn; }// modal pointer array
   inline double** GetExactSolution(){ return Qex_proj; } // continuous uniform exact solution
   inline double* GetVertexNumSol(){ return Qv; }//only numerical solution at nodes only and averaged
   inline double* GetContUniformNumSol(){ return Q_cont_sol; } // continuous uniform numerical solution
   inline double* GetContUniformExactSol(){ return Q_exact; }// continuous uniform exact solution
   inline double GetPhyTime(){ return phy_time; }

   virtual void print_cont_vertex_sol()=0;
   virtual void print_average_sol()=0;
   virtual void dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                    ,double& L1_aver_sol_,double& L2_aver_sol_
                    ,double& L1_nodal_gausspts, double& L2_nodal_gausspts)=0;
   virtual void dump_discont_sol()=0;
   virtual void dump_timeaccurate_sol()=0;
   virtual void dump_timeaccurate_errors()=0; // remember that you always have to dump the time_accurate solution first
   virtual void dump_timeaccurate_waveenergy(const double& in_E_ex_,
                                             const double& in_E_,
                                             const double& in_GG_ex_,
                                             const double& in_GG_)=0;
protected:

   virtual void UpdateResidOneCell(const int& cellid, double* q_
                           , double* resid_)=0;

   virtual double eval_init_sol(const double& xx)=0;
   virtual double eval_exact_sol(double& xx)=0;
   virtual double eval_basis_poly(const double& xi_, const int& basis_k_)=0;
   virtual double eval_basis_norm_squared(const int& basis_k_)=0;
   virtual void setup_basis_interpolation_matrices()=0;

   virtual double evalSolution(const double *q_, const double& xi_pt)=0;
   virtual double evalSolution(const double *q_, const int& position_)=0;
   virtual double initSol_legendre_proj(const int& eID, const int &basis_id,
                                const GaussQuad & quad_)=0;
   virtual double ExactSol_legendre_proj(const int &eID,
                                 const int &basis_k,
                                  const GaussQuad &quad_)=0;

   virtual void CalcTimeStep()=0;              // calculate time step and CFL
   virtual void ComputeExactSolShift()=0;      // compute shift each time step
   virtual void Compute_exact_vertex_sol()=0;  //exact continuous solution
   virtual void Compute_projected_exact_sol()=0; //projected exact solution

   virtual void compute_uniform_cont_sol()=0; // compute numerical cont solution
   virtual double compute_totalVariation()=0; // compute TV if needed

   virtual double TimeAccurateExactSol_legendre_proj(const int &eID,
                                                const int &basis_k,
                                                const GaussQuad &quad_)=0;
   virtual void Compute_TimeAccurate_exact_sol()=0;

   //virtual double Compute_localEnergyOneCell(const int& eID_)=0;
   //virtual double Compute_localEnergyOneCell(const int& eID_)=0;

protected:

   GridData *grid_=nullptr;
   SimData *simdata_=nullptr;

   int Ndof = 1;

   double **Qn=nullptr;      // Nelem * Ndof long
   double *Q_exact=nullptr;  // Nfaces long
   double **Qex_proj=nullptr; // projected exact solution , Nelem long
   double *Qv=nullptr;       // Nfaces long

   double *flux_com=nullptr;  // common interface flux, Nfaces long

   double max_eigen_advec=0.0;  // maximum eigenvalue for adevction
   double phy_time=0.0;         // physical time
   double time_step=1e-5;       // time step dt
   double last_time_step=1e-5;  // last time step dt_last
   double CFL=1.0;              // CFL number

   double exact_sol_shift=0.;  // at or adt, since exact sol have f(x-at).
   double wave_length_=0.;     // wave length, Lambda
   double wave_speed_=0.0;      // wave speed (a), u_t + a u_x =0
   double init_wave_E_=0.0;
   double wave_energy_=0.0;
   double Peclet_no=1.0;

   int Nquad_=1; // Gauss Quadrature rules

   GaussQuad quad_;

   double **Lk=nullptr;   // Legendre polynomials at (-1) and (1)
   double *Lk_norm_squar=nullptr; // norm of squared of Lk's

   double *Q_cont_sol=nullptr;  // continuous numerical solution with uniform spacing

};

#endif
