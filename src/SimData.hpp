#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"general_tools.h"
#include"global_var.h"


struct SimData {

    std::string eqn_set;
    int poly_order_=0;    // FD Scheme order
    int RK_order_=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    int print_freq_=10;
    double penalty_param_=1.0;
    std::string diffus_scheme_type_;

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
    std::string case_no_;  // case no for burgers decay turb

    double a_wave_=2;    // wave speed
    int wave_form_ = 0;  // 0: sine wave, 1: Gaussian wave
    double Gaussian_exponent_ = -40; // u(x) = exp(-38.6 *x^2)
    double wave_freq_= 2.0;
    double thermal_diffus=1.0;

    int Nelem_ = 1;  // no. of elements in the grid
    int Npplot = 1;  // no. of equally spaced points per element for plotting
    int N_exact_plot_pts=100;
    int N_uniform_pts_per_elem_ = 100; // no. of global equally spaced points for plotting

    double x0_=0.0;
    double xf_=1.0;

    int uniform_=1;  // 0: for nonuniform mesh elements

    int refine_level_=0; // 0: no refinement

    double upwind_param_=1.0;

    // Burger's Tubulence Parameters:
    std::string turb_prob_type_;
    int max_wave_no_ = 1024;
    double max_energy_wave_no_ = 10.0;
    int* k_wave_no_ =nullptr;
    double* epsi_phase_=nullptr;
    double* energy_spect_=nullptr;
    int spectrum_restart_flag = 0;
    double data_print_time_=0.01;


    void Parse(const std::string &fname);

    void setup_output_directory();
    void dump_python_inputfile();

    void prepare_dump_burgers_turb_param();

    void Reset();
};


#endif
