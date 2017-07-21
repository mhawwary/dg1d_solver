#include"SimData.hpp"

void SimData::Parse(const std::string &fname){

    GetPot gp_input(fname.c_str());

    Nelem_ = gp_input("Case/num_elements",1);
    x0_ = gp_input("Case/x_begin",0.0);
    xf_ = gp_input("Case/x_end",1.0);
    uniform_ = gp_input("Case/uniform_grid",1);
    refine_level_ = gp_input("Case/refinement_level",0);
    Npplot = gp_input("Case/Npoints_plot",2);
    N_exact_plot_pts = gp_input("Case/N_exact_plot_pts",100);

    print_freq_=gp_input("Simulation/print_freq",0);
    restart_iter_ = gp_input("Simulation/restart_iter",0);
    restart_flag = gp_input("Simulation/restart_flag",0);
    Sim_mode = gp_input("Simulation/mode","normal");

    a_wave_ = gp_input("wave/wave_speed",1.0);
    wave_form_ = gp_input("wave/wave_form",0);
    Gaussian_exponent_ = gp_input("wave/Gaussian_exponent",-50.0);
    wave_freq_ = gp_input("wave/wave_frequency",2.0);

    poly_order_=gp_input("space_solver/polynomial_order",1);
    upwind_param_=gp_input("space_solver/upwind_param",1.0);
    eqn_set = gp_input("space_solver/eqn_set","Advection");
    thermal_diffus
            = gp_input("space_solver/heat_eqn/thermal_diffusivity",1.0);
    penalty_param_
            = gp_input("space_solver/heat_eqn/penalty_param",1.0);
    diffus_scheme_type_ =
            gp_input("space_solver/heat_eqn/diffus_scheme_type","SIP");

    calc_dt_flag = gp_input("time_solver/calculate_dt_flag",1);
    CFL_ = gp_input("time_solver/CFL_no",1.0);
    dt_ = gp_input("time_solver/dt",1e-9);
    t_init_ = gp_input("time_solver/initial_time",0.0);
    t_end_ = gp_input("time_solver/final_time",1.0);
    maxIter_ = gp_input("time_solver/maximum_iteration",1e9);
    end_of_sim_flag_ = gp_input("time_solver/end_of_simulation_flag",0);
    RK_order_=gp_input("time_solver/explicit/RK_order",1);
    Nperiods = gp_input("time_solver/no_of_periods",1.0);

}

void SimData::setup_output_directory(){

    // Setting up some directories:
    //---------------------------------
    allocator<char> allchar; // default allocator for char

    case_postproc_dir =new char[300];

    char *case_dir=nullptr,*case_title=nullptr;
    case_dir=new char[70];
    case_title=new char[20];

    if(wave_form_==0) sprintf(case_title,"sine_wave");
    else if(wave_form_==1) sprintf(case_title,"Gaussian_wave");

    if(Sim_mode=="normal")
        sprintf(case_dir,"DGp%d_RK%d",poly_order_,RK_order_);
    else if(Sim_mode=="test")
        sprintf(case_dir,"DGp%d_RK%d_test",poly_order_,RK_order_);
    else if(Sim_mode=="error_analysis_CFL"
            || Sim_mode=="error_analysis_dt"
            || Sim_mode=="error_analysis_Beta")
        sprintf(case_dir,"DGp%d_RK%d_error_analysis",poly_order_,RK_order_);
    else FatalError_exit("Simulation mode is not defined");

    char *current_working_dir=allchar.allocate(1500);
    getcwd(current_working_dir,1500);

    char *main_dir=nullptr; main_dir = new char[25];

    if(eqn_set=="Advection"){
        sprintf(main_dir,"./Results");
    }else if (eqn_set=="Diffusion"){
        sprintf(main_dir ,"./Results_diffus");
    }else{
        FatalError_exit("Equation set when specifying output directory");
    }

    mkdir(main_dir,0777);

    chdir(main_dir);

    mkdir(case_title,0777);

    chdir(case_title);

    mkdir(case_dir,0777);

    case_postproc_dir = new char[200];

//    sprintf(case_postproc_dir,"./Results/%s/%s/",case_title,case_dir);
    sprintf(case_postproc_dir,"%s/%s/%s/",main_dir,case_title,case_dir);

    chdir(case_dir);

    mkdir("./aver",0777);
    mkdir("./nodal",0777);
    mkdir("./errors",0777);

    chdir(current_working_dir);

    cout<<"\n--> Currnet working directory: "<<current_working_dir<<endl;
    cout<<"--> Post processing directory: "<<case_postproc_dir<<endl;

    emptyarray(case_dir);
    emptyarray(main_dir);
    emptyarray(case_title);

    return;
}

void SimData::dump_python_inputfile(){

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./input/python_input.in");

    FILE* python_out = fopen(fname,"w");

    fprintf(python_out,"dir:%s\n",case_postproc_dir);
    fprintf(python_out,"errors:%s\n",(char*)"errors/errors");
    fprintf(python_out,"aver:%s\n",(char*)"aver/u_aver");
    fprintf(python_out,"cont_exact:%s\n",(char*)"nodal/u_cont_exact");
    fprintf(python_out,"cont_num:%s\n",(char*)"nodal/u_cont");
    fprintf(python_out,"discont:%s\n",(char*)"nodal/u_disc");
    fprintf(python_out,"discont_exact:%s\n",(char*)"nodal/u_disc_exact");
    fprintf(python_out,"DGp:%d\n",poly_order_);
    fprintf(python_out,"RK:%d\n",RK_order_);
    fprintf(python_out,"Nelem:%d\n",Nelem_);
    fprintf(python_out,"dt:%1.3e\n",dt_);
    fprintf(python_out,"CFL:%1.3f\n",CFL_);

    if(eqn_set=="Advection"){
        fprintf(python_out,"Beta:%1.2f\n",upwind_param_);
        fprintf(python_out,"Eqn_set:%s\n",(char*)"Advection");
    }else if(eqn_set=="Diffusion"){
        fprintf(python_out,"Epsilon:%1.2f\n",penalty_param_);
        fprintf(python_out,"Eqn_set:%s\n",(char*)"Diffusion");
        fprintf(python_out,"Diffusion_scheme:%s\n",diffus_scheme_type_.c_str());
    }else{
        FatalError_exit("Eqn_set not defined in python_dump");
    }

    fprintf(python_out,"T:%1.3f\n",Nperiods);

    if(Sim_mode=="error_analysis_dt")
        fprintf(python_out,"mode:%s","dt_const");
    else if(Sim_mode=="error_analysis_CFL")
        fprintf(python_out,"mode:%s","CFL_const");
    else if(Sim_mode=="error_analysis_Beta")
        fprintf(python_out,"mode:%s","Beta_const");
    else if(Sim_mode=="error_analysis_Epsilon")
        fprintf(python_out,"mode:%s","Epsilon_const");
    else
        fprintf(python_out,"mode:%s",Sim_mode.c_str());

    fclose(python_out);
    emptyarray(fname);

    return;
}

