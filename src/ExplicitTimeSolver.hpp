#include"DGSolver.hpp"
#include"SimData.hpp"

class ExplicitTimeSolver{

public:
    ExplicitTimeSolver(void);
    ~ExplicitTimeSolver(void);
    void setupTimeSolver(DGSolver* dg_solver_, SimData* simdata_);
    void SolveOneStep(double **qn_);
    int GetIter();

public:
    DGSolver *space_solver=nullptr;
    SimData  *simdata=nullptr;

protected:

    void FwdEuler(double **q_);
    void SSPRK22(double **q_);
    void SSPRK33(double **q_);

    void CopyOldSol(double **q_t_, double **qn_);

    void UpdateIter();

    void Reset_time_solver();

protected:
    double **resid=nullptr;
    double **q_temp=nullptr;

    int Ndof=1;

    int IterNo=0;

    double dt_=0.0;

};
