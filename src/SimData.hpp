#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"general_tools.h"


struct SimData {

    int poly_order_=0;    // FD Scheme order
    int RK_order_=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    int print_freq_=10;

    char *case_postproc_dir=nullptr;

    int calc_dt_flag=1; // 1: specify CFL and calc dt, 0: specify dt and calc CFL
    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0.0;  // initial time
    double t_end_ =1e20;  // end time
    double maxIter_ = 1e10; // maximum number of iterations
    double CFL_    = 1.0;   // CFL no.
    double Nperiods = 1.0; // no. of periods for simulation

    int restart_flag=0;  //0: start a new simulation; 1: restart a previous simulation
    int restart_iter_=0;
    int end_of_sim_flag_=0;  // 1: use max_iteration as a stopping criteria if not converged or diverged
    std::string Sim_mode;

    double a_wave_=2;    // wave speed
    int wave_form_ = 0;  // 0: sine wave, 1: Gaussian wave
    double Gaussian_exponent_ = -40; // u(x) = exp(-38.6 *x^2)

    int Nelem_ = 1;  // no. of elements in the grid
    int Npplot = 1;  // no. of equally spaced points per element for plotting

    double x0_=0.0;
    double xf_=1.0;

    int uniform_=1;  // 0: for nonuniform mesh elements

    int refine_level_=0; // 0: no refinement

    double upwind_param_=1.0;

    void Parse(const std::string &fname);

    void setup_output_directory();
    void dump_python_inputfile();
};


#endif
