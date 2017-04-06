#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    Nelem_ = gp_input("Case/num_elements",1);

    x0_ = gp_input("Case/x_begin",0.0);

    xf_ = gp_input("Case/x_end",1.0);

    uniform_ = gp_input("Case/uniform_grid",1);

    refine_level_ = gp_input("Case/refinement_level",0);


    print_freq_=gp_input("Simulation/print_freq",0);


    a_wave_ = gp_input("wave/wave_speed",1.0);


    poly_order_=gp_input("space_solver/polynomial_order",1);
    upwind_param_=gp_input("space_solver/upwind_param",1);

    calc_dt_flag = gp_input("time_solver/calculate_dt_flag",1);

    CFL_ = gp_input("time_solver/CFL_no",1.0);

    dt_ = gp_input("time_solver/dt",1e-9);

    t_init_ = gp_input("time_solver/initial_time",0.0);

    t_end_ = gp_input("time_solver/final_time",1.0);

    maxIter_ = gp_input("time_solver/maximum_iteration",1e9);


    RK_order_=gp_input("time_solver/explicit/RK_order",1);

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

