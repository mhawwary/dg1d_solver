#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"general_tools.h"


struct SimData {

    int restart_flag=0;  //0: start a new simulation; 1: restart a previous simulation
    int poly_order_=0;    // FD Scheme order
    int RK_order_=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    int print_freq_=10;

    int calc_dt_flag=1; // 1: specify CFL and calc dt, 0: specify dt and calc CFL
    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0.0;  // initial time
    double t_end_ =1e20;  // end time
    double maxIter_ = 1e10; // maximum number of iterations
    double CFL_    = 1.0;   // CFL no.
    double a_wave_=2;

    int  Nelem_ = 1;  // no. of elements in the grid

    double x0_=0.0;
    double xf_=1.0;

    int uniform_=1;  // 0: for nonuniform mesh elements

    int refine_level_=0; // 0: no refinement

    int upwind_param_=1;

    void Parse(const std::string &fname);

    void print_data();

};


#endif
