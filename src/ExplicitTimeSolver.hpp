#include"DGSolver.hpp"
#include"SimData.hpp"

class ExplicitTimeSolver{

public:
    ExplicitTimeSolver(void);
    ~ExplicitTimeSolver(void);
    void setupTimeSolver(DGSolver* dg_solver_, SimData* simdata_);
    void SolveOneStep(double **qn_);
    void ComputeInitialResid(double **qn_);

    int GetIter(){

        return IterNo;
    }

    void Set_time_step(const double& dt_set){
        dt_ = dt_set;
        return;
    }

    void Reset_iter(const double& iter_reset){
        IterNo = iter_reset;
        return;
    }

public:
    DGSolver *space_solver=nullptr;
    SimData  *simdata=nullptr;

protected:

    void FwdEuler(double **q_);
    void SSPRK22(double **q_);
    void SSPRK33(double **q_);
    void classicRK4(double **q_);
    //void RK54(double **q_);

    void CopyOldSol(double **q_t_, double **qn_);
    void CopyOldResid(double **resid_t_, double **old_resid_);

    void UpdateIter(){

        IterNo++;

        return;
    }

    void Reset_time_solver();

protected:
    double **resid=nullptr;
    double **q_temp=nullptr;
    double **resid_temp0=nullptr;
    double **resid_temp1=nullptr;
    double **resid_temp2=nullptr;

    int Ndof=1;

    int IterNo=0;

    double dt_=0.0;

};
