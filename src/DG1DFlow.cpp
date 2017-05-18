#include "general_tools.h"
#include "SimData.hpp"
#include "GirdData.h"
#include "DGSolver.hpp"
#include "ExplicitTimeSolver.hpp"


SimData simdata_;
GridData meshdata_;
DGSolver dg_solver_;
ExplicitTimeSolver time_solver_;

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

    return 0;
}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata_.Parse(argv[argc-1]);
        simdata_.setup_output_directory();
        simdata_.dump_python_inputfile();
    }

    meshdata_.set_grid_param(simdata_);

    meshdata_.generate_grid();

    dg_solver_.setup_solver(meshdata_,simdata_);

    dg_solver_.InitSol();

    return;
}

void RunSim(){

    time_solver_.setupTimeSolver(&dg_solver_,&simdata_);

    double gtime = dg_solver_.GetPhyTime();

    double dt_= dg_solver_.GetTimeStep();

    time_solver_.ComputeInitialResid(dg_solver_.GetNumSolution());

    time_solver_.SolveOneStep(dg_solver_.GetNumSolution());

    time_solver_.space_solver->UpdatePhyTime(dt_);

    gtime=dg_solver_.GetPhyTime();

    while ( gtime < (simdata_.t_end_-0.5*dt_) ){

            time_solver_.SolveOneStep(dg_solver_.GetNumSolution());

            time_solver_.space_solver->UpdatePhyTime(dt_);

            gtime=dg_solver_.GetPhyTime();
    }

    return;
}

void PostProcess(){

    double aver_L2_error=0.0,proj_sol_L2_error=0.0;

    dg_solver_.Compute_vertex_sol();
    proj_sol_L2_error = dg_solver_.Compute_projected_sol_error();
    aver_L2_error= dg_solver_.ComputeAverageError();

    dg_solver_.print_cont_vertex_sol();
    dg_solver_.print_average_sol();
    dg_solver_.dump_discont_sol();
    dg_solver_.dump_errors(proj_sol_L2_error,aver_L2_error);

    printf("\nFinal Iteration number is: %d\n",time_solver_.GetIter());
    printf("Final time is: %1.2f\n",dg_solver_.GetPhyTime());
    printf("proj_sol_L2_Error: %e  , Aver_L2_error: %e\n\n"
           ,proj_sol_L2_error,aver_L2_error);


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












