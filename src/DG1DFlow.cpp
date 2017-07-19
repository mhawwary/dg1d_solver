#include "general_tools.h"
#include "SimData.hpp"
#include "GirdData.h"
#include "DGSolver.hpp"
#include "DGSolverAdvec.hpp"
#include "DGSolverDiffus.hpp"
#include "ExplicitTimeSolver.hpp"


SimData simdata_;
GridData meshdata_;
DGSolver *dg_solver_;
ExplicitTimeSolver *time_solver_;

void InitSim(const int& argc, char** argv);
void RunSim();
void PostProcess();
void logo();

int main(int argc, char** argv){

    if (argc < 2) {

        FatalErrorST("ERROR: No inputs are specified ... ");

        return(0);
    }

    logo();

    InitSim(argc, argv);

    RunSim();

    PostProcess();

    emptypointer(dg_solver_);
    emptypointer(time_solver_);

    return 0;
}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata_.Parse(argv[argc-1]);
        simdata_.setup_output_directory();
    }

    meshdata_.set_grid_param(simdata_);

    meshdata_.generate_grid();

    // Allocating Solvers:
    if(simdata_.eqn_set=="Advection")
        dg_solver_ = new DGSolverAdvec;
    else if(simdata_.eqn_set=="Diffusion")
        dg_solver_ = new DGSolverAdvec;

    time_solver_ = new ExplicitTimeSolver;

    dg_solver_->setup_solver(meshdata_,simdata_);

    dg_solver_->InitSol();

    time_solver_->setupTimeSolver(dg_solver_,&simdata_);

    simdata_.dump_python_inputfile();

    return;
}

void RunSim(){

    double gtime = dg_solver_->GetPhyTime();

    double dt_= dg_solver_->GetTimeStep();

    time_solver_->ComputeInitialResid(dg_solver_->GetNumSolution());

    time_solver_->SolveOneStep(dg_solver_->GetNumSolution());

    time_solver_->space_solver->UpdatePhyTime(dt_);

    gtime=dg_solver_->GetPhyTime();

    while ( gtime < simdata_.t_end_- 1.05*dt_ ){

        time_solver_->SolveOneStep(dg_solver_->GetNumSolution());

        time_solver_->space_solver->UpdatePhyTime(dt_);

        gtime=dg_solver_->GetPhyTime();
    }

    // Last iteration:

    time_solver_->SolveOneStep(dg_solver_->GetNumSolution());

    time_solver_->space_solver->UpdatePhyTime(dg_solver_->GetLastTimeStep());

    gtime=dg_solver_->GetPhyTime();

    return;
}

void PostProcess(){

    double L1_aversol_=0.0,L1_projsol_=0.0
            ,L2_aversol_=0.0,L2_projsol_=0.0
            , L1_nodal_gausspts, L2_nodal_gausspts;

    dg_solver_->Compute_vertex_sol();
    L1_projsol_ = dg_solver_->L1_error_projected_sol();
    L2_projsol_ = dg_solver_->L2_error_projected_sol();
    L1_aversol_ = dg_solver_->L1_error_average_sol();
    L2_aversol_ = dg_solver_->L2_error_average_sol();
    L1_nodal_gausspts = dg_solver_->L1_error_nodal_gausspts();
    L2_nodal_gausspts = dg_solver_->L2_error_nodal_gausspts();

    dg_solver_->print_cont_vertex_sol();
    dg_solver_->print_average_sol();
    dg_solver_->dump_discont_sol();
    dg_solver_->dump_errors(L1_projsol_,L2_projsol_
                           ,L1_aversol_, L2_aversol_
                           ,L1_nodal_gausspts, L2_nodal_gausspts);

    printf("\nFinal Iteration number is: %d\n",time_solver_->GetIter());
    printf("Final time is: %1.5f\n",dg_solver_->GetPhyTime());
    printf("L1_error,  proj_sol: %e   ,  aver_sol: %e , nodal_sol: %e\n"
           ,L1_projsol_, L1_aversol_, L1_nodal_gausspts);
    printf("L2_error,  proj_sol: %e   ,   aver_sol: %e , nodal_sol: %e\n"
           ,L2_projsol_, L2_aversol_, L2_nodal_gausspts);

    return;
}

void logo(){

    cout<<"_________________________________________________________________________________________"<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                 "<<"  Welcome to the Discontinuous Galerkin solver  "<<"                "<<endl;
    cout<<"                  "<<"   for 1D wave and scalar conservation laws    "<<"                "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"       Author:               Mohammad Alhawwary, PhD. Student                            "<<endl;
    cout<<"  Affiliation:   Aerospace Engineering Department, University of Kansas, USA             "<< endl;
    cout<<"_______________________________________04/05/2017________________________________________"<<endl;

    return;
}












