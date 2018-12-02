#include "general_tools.h"
#include "SimData.hpp"
#include "GridData.h"
#include "DGSolver.hpp"
#include "DGSolverAdvec.hpp"
#include "DGSolverDiffus.hpp"
#include "DGSolverAdvecDiffus.hpp"
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

    clock_t t_start=clock();

    InitSim(argc, argv);

    RunSim();

    clock_t t_end=clock();

    _print_log("\n\nEnd simulation");
    cout << "Elapsed Time: " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC
         << " seconds\n" <<endl;

    emptypointer(dg_solver_);
    emptypointer(time_solver_);

    return 0;
}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata_.Parse(argv[argc-1]);
        simdata_.setup_output_directory();

        if(simdata_.wave_form_==3) //Burgers Turbulence
            simdata_.prepare_dump_burgers_turb_param();
    }

    _print_log("Init grid");
    meshdata_.set_grid_param(simdata_);
    meshdata_.generate_grid();

    // Allocating Solvers:
    _print_log("Init space solver")
    if(simdata_.eqn_set=="Advection")
        dg_solver_ = new DGSolverAdvec;
    else if(simdata_.eqn_set=="Diffusion")
        dg_solver_ = new DGSolverDiffus;
    else if(simdata_.eqn_set=="Advection_Diffusion")
        dg_solver_ = new DGSolverAdvecDiffus;
    else
        _notImplemented("Equation set");

    dg_solver_->setup_solver(meshdata_,simdata_);
    dg_solver_->InitSol();

    _print_log("Init time solver")
    time_solver_ = new ExplicitTimeSolver;
    time_solver_->setupTimeSolver(dg_solver_,&simdata_);
    simdata_.dump_python_inputfile();

    return;
}

void RunSim(){

    int n_iter_print;
    int local_iter=0;
    double gtime = dg_solver_->GetPhyTime();
    double dt_= dg_solver_->GetTimeStep();
    double dt_last_print=0.0;
    double temp_tol=1e-8;

    //======================================================================
    //             Preparing simulation control variables
    //======================================================================
    if(simdata_.unsteady_data_print_flag_==0){    // use iter_print to print
        n_iter_print = simdata_.unsteady_data_print_iter_;
        if(n_iter_print<=1)
            FatalError_exit("Warning: iter to print is very small <=1 ");
        dt_last_print = dt_;
        n_iter_print--;

    }else if(simdata_.unsteady_data_print_flag_==1){  // use time point to print
        if(simdata_.unsteady_data_print_time_ < (dt_ - temp_tol))
            FatalError_exit("Warning:  time to print is less than dt");

        n_iter_print= (int) round( simdata_.unsteady_data_print_time_/ dt_) ;

        if((n_iter_print*dt_) > (simdata_.unsteady_data_print_time_-temp_tol) ){
            n_iter_print--;
            dt_last_print = simdata_.unsteady_data_print_time_ - (n_iter_print * dt_);
        }else if((n_iter_print*dt_) < (simdata_.unsteady_data_print_time_+temp_tol) ){
            dt_last_print = simdata_.unsteady_data_print_time_ - (n_iter_print*dt_);
        }

        if(n_iter_print<=1)
            FatalError_exit("Warning: iter to print is very small <=1 ");

    // print using the specified iter_print without dt changing except the last one
    }else if(simdata_.unsteady_data_print_flag_==2){
        n_iter_print=simdata_.unsteady_data_print_iter_;
        dt_last_print = dg_solver_->GetLastTimeStep();
    }else{
        FatalError_exit("unsteady data print flag error");
    }

    printf("N_iter to print unsteady data: %d, dt_last: %1.5e\n"
           ,n_iter_print, dt_last_print);

    printf("Iter No:%d, time: %1.5f",time_solver_->GetIter()
           ,dg_solver_->GetPhyTime());
    // Dump initial data:
    dg_solver_->dump_timeaccurate_sol();
    dg_solver_->dump_timeaccurate_errors();
    //===========================
    // Solve First Iteration
    //===========================
    _print_log("\nComputing initial residual & solve one step");
    time_solver_->ComputeInitialResid(dg_solver_->GetNumSolution());
    time_solver_->SolveOneStep(dg_solver_->GetNumSolution());
    time_solver_->space_solver->UpdatePhyTime(dt_);
    gtime=dg_solver_->GetPhyTime();
    local_iter++;

    if(n_iter_print==1){
        printf("\nIter No:%d, time: %f",time_solver_->GetIter(),gtime);
        dg_solver_->dump_timeaccurate_sol();
        dg_solver_->dump_timeaccurate_errors();
        local_iter=0;
    }

    //======================================================================
    //                        Main Solution Loop
    //======================================================================
    _print_log("\nBegin main solution loop (time advance)");
    if(simdata_.unsteady_data_print_flag_==0
            || simdata_.unsteady_data_print_flag_==1){
        while ( gtime < (simdata_.t_end_-(1+1e-5)*(dt_+1e-10))){
            time_solver_->SolveOneStep(dg_solver_->GetNumSolution());
            time_solver_->space_solver->UpdatePhyTime(dt_);
            gtime=dg_solver_->GetPhyTime();
            local_iter++;

            if(local_iter%n_iter_print==0){
                time_solver_->Set_time_step(dt_last_print);
                time_solver_->SolveOneStep(dg_solver_->GetNumSolution());
                time_solver_->space_solver->UpdatePhyTime(dt_last_print);
                gtime=dg_solver_->GetPhyTime();
                printf("\nIter No:%d, time: %1.5f",time_solver_->GetIter(),gtime);
                dg_solver_->dump_timeaccurate_sol();
                dg_solver_->dump_timeaccurate_errors();
                time_solver_->Set_time_step(dt_);
                //time_solver_->Reset_iter(time_solver_->GetIter()-1);
                local_iter=0;
            }
        }

    }else if(simdata_.unsteady_data_print_flag_==2){
        while ( fabs(gtime - simdata_.t_end_) > (dt_+temp_tol) ){

            time_solver_->SolveOneStep(dg_solver_->GetNumSolution());
            time_solver_->space_solver->UpdatePhyTime(dt_);
            gtime=dg_solver_->GetPhyTime();
            local_iter++;

            if(local_iter%n_iter_print==0){
                printf("\nIter No:%d, time: %1.5f",time_solver_->GetIter(),gtime);
                dg_solver_->dump_timeaccurate_sol();
                dg_solver_->dump_timeaccurate_errors();
                local_iter=0;
            }
        }

        // Last iteration:
        dt_last_print = dg_solver_->GetLastTimeStep();
        time_solver_->Set_time_step(dt_last_print);
        time_solver_->SolveOneStep(dg_solver_->GetNumSolution());
        if(dt_last_print>=temp_tol)
            time_solver_->space_solver->UpdatePhyTime(dt_last_print);
        else
            time_solver_->space_solver->UpdatePhyTime(dt_);

        gtime=dg_solver_->GetPhyTime();
        printf("\nIter No:%d, time: %1.5f",time_solver_->GetIter(),gtime);
        dg_solver_->dump_timeaccurate_sol();
        dg_solver_->dump_timeaccurate_errors();
    }

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
    cout<<"                "<<"  Welcome to the Discontinuous Galerkin solver   "<<"                "<<endl;
    cout<<"                "<<"  for 1D wave/heat and scalar conservation laws  "<<"                "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"       Author:               Mohammad Alhawwary, PhD. Student                            "<<endl;
    cout<<"       Affiliation:  Aerospace Engineering Department, University of Kansas, USA         "<< endl;
    cout<<"_______________________________________04/05/2017________________________________________"<<endl;

    return;
}










