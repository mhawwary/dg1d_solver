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

int main(int argc, char** argv){

    if (argc < 2) {

        FatalErrorST("ERROR: No inputs are specified ... ");

        return(0);
    }

    InitSim(argc, argv);

    RunSim();

    _print("back to the main function");

    PostProcess();

    _print("finished post processing");

    return 0;
}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata_.Parse(argv[argc-1]);
        simdata_.print_data();
    }

    meshdata_.set_grid_param(simdata_);

    cout << "\n--finished setting up grid parameters\n";

    meshdata_.generate_grid();

    cout << "\n--finished generating the grid\n";

    dg_solver_.setup_solver(meshdata_,simdata_);

    dg_solver_.InitSol();

    cout << "\n--finished Initializing the Solution\n";

    return;
}

void RunSim(){

    int n=0;

    double gtime=dg_solver_.GetPhyTime();

    double dt_= dg_solver_.GetTimeStep();

    time_solver_.setupTimeSolver(&dg_solver_,&simdata_);

    _print("finished setting up time solver");

    while ( gtime <=fabs( simdata_.t_end_ - pow(10,-10)) ){

            gtime += dt_;

            time_solver_.space_solver->UpdatePhyTime(dt_);

            time_solver_.SolveOneStep(dg_solver_.GetNumSolution());

            n=time_solver_.GetIter();
        }

    _(n);
    _(gtime);
    _print("finished time loop");

    return;
}

void PostProcess(){

    meshdata_.print_grid();

    dg_solver_.Compute_vertex_sol();
    dg_solver_.Compute_exact_sol();

    dg_solver_.print_exact_sol();
    dg_solver_.print_num_vertex_sol();

    dg_solver_.print_num_average_sol();
    dg_solver_.print_exact_average_sol();

    return;
}














