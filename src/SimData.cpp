#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    Nelem_ = gp_input("Case/num_elements",1);

    x0_ = gp_input("Case/x_begin",0.0);

    xf_ = gp_input("Case/x_end",1.0);

    uniform_ = gp_input("Case/uniform_grid",1);

    refine_level_ = gp_input("Case/refinement_level",0);


    print_freq_=gp_input("Simulation/print_freq",0);
    restart_iter_ = gp_input("Simulation/restart_iter",0);
    restart_flag = gp_input("Simulation/restart_flag",0);


    a_wave_ = gp_input("wave/wave_speed",1.0);
    wave_form_ = gp_input("wave/wave_form",0);
    Gaussian_exponent_ = gp_input("wave/Gaussian_exponent",-50.0);


    poly_order_=gp_input("space_solver/polynomial_order",1);
    upwind_param_=gp_input("space_solver/upwind_param",1.0);

    calc_dt_flag = gp_input("time_solver/calculate_dt_flag",1);

    CFL_ = gp_input("time_solver/CFL_no",1.0);

    dt_ = gp_input("time_solver/dt",1e-9);

    t_init_ = gp_input("time_solver/initial_time",0.0);

    t_end_ = gp_input("time_solver/final_time",1.0);

    maxIter_ = gp_input("time_solver/maximum_iteration",1e9);

    end_of_sim_flag_ = gp_input("time_solver/end_of_simulation_flag",0);

    RK_order_=gp_input("time_solver/explicit/RK_order",1);

    Nperiods = gp_input("time_solver/no_of_periods",1);


    // Setting up some directories:
    //---------------------------------
    allocator<char> allchar; // default allocator for char

    case_postproc_dir =new char[200];

    char *case_dir=nullptr; case_dir=new char[100];

    sprintf(case_dir,"DGp%d_RK%d",poly_order_,RK_order_);

    char *current_working_dir=allchar.allocate(1500);
    getcwd(current_working_dir,1500);

    chdir("./Results");

    mkdir(case_dir,0777);

    case_postproc_dir = new char[200];

    sprintf(case_postproc_dir,"./Results/%s/",case_dir);

    chdir(current_working_dir);

    cout<<"\n--> Currnet working directory: "<<current_working_dir<<endl;
    cout<<"--> Post processing directory: "<<case_postproc_dir<<endl;

    emptyarray(case_dir);

}

void SimData::print_data(){

    // Screen Output of input and simulation parameters:
    cout <<"\n===============================================\n";
    cout << "CFL no.:  "<<CFL_<<"\tWave Speed:  "<<a_wave_<<endl;
    cout << "dt:  "<<dt_<<"\t"<< "dx:  "<<a_wave_*dt_/CFL_<<endl;
    cout << "required no. of time steps: "<<t_end_/dt_<<endl;
    cout << "Number of Elements:  "<<Nelem_<<endl;
    cout << "Polynomial Order: "<<poly_order_<<endl;
    cout << "RK_order:  "<< RK_order_ << endl <<"\n";

    return;
}

