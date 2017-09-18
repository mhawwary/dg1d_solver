#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"general_tools.h"
#include"global_var.h"


struct SimData {

    // Case parameters:
    //-----------------------------
    std::string case_title;    // parent folder for the case results
    int Nelem_ = 1;  // no. of elements in the grid
    double x0_=0.0;  // x_begin of domain
    double xf_=1.0;  // x_end of domain
    int uniform_=1;  // 0: for nonuniform mesh elements
    int refine_level_=0; // 0: no refinement
    int Npplot = 1;  // to be depricated to use N_uniform_pts_per_elem_ instead
    int N_exact_plot_pts=100;   // no. of points for exact solution plotting
    int N_uniform_pts_per_elem_ = 100; // no. of global equally spaced points for plotting

    // Simulation parameters:
    //-----------------------------
    int unsteady_data_print_flag_ = 1;   // 0: use iter , 1: use time
    int unsteady_data_print_iter_ = 1000;  // iter no. for print of unsteady data
    double unsteady_data_print_time_ = 1000; // time point for print of unsteady data
    int restart_flag=0;      //  0: start a new simulation; 1: restart a previous simulation
    int restart_iter_=0;
    std::string Sim_mode;    // normal/test/error_analysis_dt/error_analysis_CFL/CFL_const/dt_const
    std::string case_no_;    // case no for burgers decay turb

    // Wave parameters:
    //----------------------------
    double a_wave_=2;    // wave speed
    int wave_form_ = 0;  // 0: sine wave, 1: Gaussian wave
    // ./trigonometric , sine or cosine ,  u(x,0) = A * sin ( f * PI + phy ) + C :
    std::string trig_wave_type_;  // sine or cosine
    double wave_freq_= 2.0;  // f
    double wave_amp_ = 1.0;  // A
    double wave_const = 0.0; // C
    double wave_shift = 0.0; // phy
    // ./Gaussian , u(x) = A * exp( M *x^2) :
    double Gaussian_amp_ = 1.0;   // A
    double Gaussian_exponent_ = -40;  // M
    // ./Burger_turb:  E(k) = A * exp(-(k/ko)^(2)) , A = ko^(-5) / sqrt(3*PI)
    std::string turb_prob_type_;  // Decay_turb_Adams, Decay_turb_Yanan, forced_turb_Sherwin
    int max_wave_no_ = 1024;     // k_max
    double max_energy_wave_no_ = 10.0;  // ko
    int* k_wave_no_ =nullptr;           // k array
    double* epsi_phase_=nullptr;        // phy phase for u(x) based on E
    double* energy_spect_=nullptr;      // E
    int spectrum_restart_flag = 0;      // 0: compute new spectrum, 1: load spectrum
    double velocity_mean_=0.0;   // u(x) = sum(E(k) *cos(k))+ u_mean

    // Space Solver parameters:
    //-----------------------------
    std::string eqn_set;    //  Advection / Diffusion / Advection_Diffusion
    std::string eqn_type_;  // inv_burger / linear_advec / linear_diffus
    int poly_order_=0;      // polynomial order
    // ./advec_eqn:
    double upwind_param_=1.0;     // Beta for upwind-central
    // ./heat_eqn:
    double penalty_param_=1.0;    // Epsilon for relaxation
    double thermal_diffus=1.0;    // kappa
    std::string diffus_scheme_type_;   // BR2 / SIP / .....

    // Time Solver parameters:
    //--------------------------------
    int calc_dt_flag=1; // 1: specify CFL and calc dt, 0: specify dt and calc CFL
    int calc_dt_adv_diffus_flag=0;  // 0: based on advection, 1: based on diffusion, 2: based on combined advection-diffusion
    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0.0;  // initial time
    double t_end_ =1e20;  // end time
    double maxIter_ = 1e10; // maximum number of iterations
    double CFL_    = 1.0;   // CFL no.
    double Nperiods = 1.0; // no. of periods for simulation
    int end_of_sim_flag_=0;  // 1: use max_iteration as a stopping criteria if not converged or diverged
    // ./explicit:
    int RK_order_=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    //------------------------------------------------------------------

    char *case_postproc_dir=nullptr;

    // Function members:
    //--------------------------
    void Parse(const std::string &fname);

    void setup_output_directory();
    void dump_python_inputfile();

    void prepare_dump_burgers_turb_param();

    void Reset();
};


#endif
