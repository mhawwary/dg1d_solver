#include "DGSolverAdvecDiffus.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

DGSolverAdvecDiffus::~DGSolverAdvecDiffus(void){

    Reset_solver();
}

void DGSolverAdvecDiffus::setup_solver(GridData& meshdata_, SimData& osimdata_){

    register int i;

    grid_ = &meshdata_;
    simdata_ = &osimdata_;

    Ndof = simdata_->poly_order_+1;

    // Nquad  is for error integrations
    // Nquad_invFlux & Nquad_viscFlux_ are for flux projection
    switch (Ndof) {
    case 1:  // p0
        Nquad_ = 1;
        Nquad_invFlux_= 1; // it is not applicable
        Nquad_viscFlux_=1; // it is not applicable
        break;
    case 2: // p1
        Nquad_ = 2;
        Nquad_invFlux_= 2;
        Nquad_viscFlux_=1;
        break;
    case 3:  // p2
        Nquad_ = 3;
        Nquad_invFlux_= 3;
        Nquad_viscFlux_=2;
        break;
    case 4:  // p3
        Nquad_ = 4;
        Nquad_invFlux_= 5;
        Nquad_viscFlux_=3;
        break;
    case 5:  // p4
        Nquad_ = 5;
        Nquad_invFlux_= 6;
        Nquad_viscFlux_=4;
        break;
    case 6:  // p5
        Nquad_ = 6;
        Nquad_invFlux_= 8;
        Nquad_viscFlux_=5;
        break;
    default:
        break;
    }

    quad_.setup_quadrature(Nquad_);
    quad_invF_.setup_quadrature(Nquad_invFlux_);
    quad_viscF_.setup_quadrature(Nquad_viscFlux_);

    Lk = new double*[Ndof]; // Lk[k, xi]
    dLk = new double*[Ndof]; // dLk[k, xi]
    Lk_norm_squar = new double[Ndof];

    for(i=0; i<Ndof; i++){
        Lk[i] = new double[2];
        dLk[i] = new double[2];
    }

    Qn    =  new double* [grid_->Nelem];
    Qex_proj = new double*[grid_->Nelem];

    for(i=0; i<grid_->Nelem; i++){
        Qn[i]       = new double[Ndof];
        Qex_proj[i] = new double[Ndof];
    }

    Q_exact = new double[grid_->N_exact_ppts];
    Qv = new double[grid_->Nfaces];

    flux_com = new double[grid_->Nfaces];
    viscflux_com = new double[grid_->Nfaces];
    u_sol_jump = new double[grid_->Nfaces+2]; // we added 2 to account for periodic B.C. and non compactness of BR1
    Q_cont_sol = new double[grid_->N_uniform_pts];

    eta_penalty = simdata_->penalty_param_; // penalty param
    C_lift = pow(Ndof,2)/2.0; // C lifting term multiplier for LDG, BR1, BR2/SIP
    C1_lift = pow(-1,Ndof-1)*Ndof/4.0; // C1 lifting term multiplier for BR1
    if(simdata_->diffus_scheme_type_=="SIP"
            || simdata_->diffus_scheme_type_=="BR2"){
        C_lift = eta_penalty * C_lift;
        C1_lift=0.0;
    }else if(simdata_->diffus_scheme_type_=="BR1"){
        C_lift = (1+eta_penalty)*C_lift;
    }else if(simdata_->diffus_scheme_type_=="LDG"){
        C_lift = 2 * C_lift + eta_penalty;
        C1_lift=0.0;
    }

    SetPhyTime(simdata_->t_init_);
    setup_basis_interpolation_matrices();

    //Wave data:
    wave_length_ = grid_->xf - grid_->x0 ;
    wave_speed_ = simdata_->a_wave_;
    ComputeExactSolShift();

    // Computing exact solutions at time 0
    Compute_projected_exact_sol();
    Compute_exact_vertex_sol();

    return;
}

void DGSolverAdvecDiffus::setup_basis_interpolation_matrices(){

    int k;

    for(k=0; k<Ndof; k++){
        Lk[k][0] = eval_basis_poly(-1,k);
        Lk[k][1] = eval_basis_poly(1,k);
        dLk[k][0] = eval_basis_poly_derivative(-1,k);
        dLk[k][1] = eval_basis_poly_derivative(1,k);
        Lk_norm_squar[k] = eval_basis_norm_squared(k);
    }

    return;
}

void DGSolverAdvecDiffus::Reset_solver(){

    emptyarray(grid_->Nelem,Qn);
    emptyarray(Q_exact);
    emptyarray(flux_com);
    emptyarray(viscflux_com);
    emptyarray(u_sol_jump);
    emptyarray(Q_cont_sol);
    emptyarray(Qv);
    emptyarray(grid_->Nelem,Qex_proj);
    emptyarray(Ndof,Lk);
    emptyarray(Ndof,dLk);
    emptyarray(Lk_norm_squar);

    grid_->Reset_();
    simdata_->Reset();

    quad_.Reset_quad();
    quad_invF_.Reset_quad();
    quad_viscF_.Reset_quad();

    return;
}


// Solver functions
//-------------------------------------------

void DGSolverAdvecDiffus::CalcTimeStep(){

    compute_uniform_cont_sol();
    double TV_ = compute_totalVariation();

    double dx = grid_->dx;
    double dx2 = pow(dx,2);
    double radius_advec_=0.0, radius_diffus_ =0.0;

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;
    radius_advec_ =  max_eigen_advec / dx  ;
    radius_diffus_ = simdata_->thermal_diffus / dx2 ;

    if(simdata_->calc_dt_flag==1){   // use CFL as input
        CFL = simdata_->CFL_;
        if(simdata_->calc_dt_adv_diffus_flag==0)   // based on advection effect only
            time_step = CFL / radius_advec_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==1)  // based on diffusion effect only
            time_step = CFL /  radius_diffus_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==2)  // based on combined advection and diffusion effects
            time_step = CFL / ( radius_advec_ + radius_diffus_ );
        else
            FatalError_exit("Wrong Calc dt adv diffus flag");

        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){   // use dt as input

        time_step = simdata_->dt_;
        last_time_step = time_step;

        if(simdata_->calc_dt_adv_diffus_flag==0)
            CFL = time_step * radius_advec_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==1)
            CFL = time_step * radius_diffus_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==2)
            CFL = time_step * ( radius_advec_ + radius_diffus_ );
        else
            FatalError_exit("Wrong Calc dt adv diffus flag");

        simdata_->CFL_ = CFL;

    }else{

        FatalError_exit("Wrong Calc_dt_flag");
    }

    // Determining end of simulation parameters:
    //----------------------------------------------------
    double temp_tol=1e-8;
    if(simdata_->end_of_sim_flag_==0){  // use Nperiods as stopping criteria

        simdata_->t_end_ = simdata_->Nperiods * T_period;
        simdata_->maxIter_ = (int) floor(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step)
                > (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step)
                 < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==1){ // use final time as stopping criteria

        simdata_->Nperiods = simdata_->t_end_/T_period;
        simdata_->maxIter_ = (int) floor(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step)
                > (simdata_->t_end_-temp_tol) ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step)
                 < (simdata_->t_end_+temp_tol) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

        if(last_time_step<=1e-10){
            last_time_step=time_step;
            simdata_->maxIter_--;
        }

    }else if(simdata_->end_of_sim_flag_==2){  // use Max Iter as stopping criteria

        simdata_->t_end_ = simdata_->maxIter_ * time_step;
        simdata_->Nperiods = simdata_->t_end_/T_period;

    }else{
        FatalError_exit("Wrong end_of_simulation_flag");
    }

    // Screen Output of input and simulation parameters:
    cout <<"\n===============================================\n";
    cout << "max eigenvalue : "<<max_eigen_advec<<endl;
    cout << "TotalVariation : "<<TV_<<endl;
    cout << "ThermalDiffusiv: "<< simdata_->thermal_diffus << endl;
    cout << "CFL no.        : "<<CFL<<endl;
    cout << "time step, dt  : "<<time_step<<endl;
    cout << "last_time_step : "<<last_time_step<<endl;
    cout << "input Nperiods : "<<simdata_->Nperiods<<endl;
    cout << "new   Nperiods : "<<simdata_->t_end_/T_period<<endl;
    cout << "exact_sol_shift: "<<exact_sol_shift<<endl;
    cout << "T_period       : "<<T_period<<endl;
    printf("actual_end_time:%1.5f",simdata_->t_end_);
    cout <<"\nMax_iter: "<<simdata_->maxIter_<<endl;

    cout << "\nNumber of Elements: "<< grid_->Nelem<<"  dx:  "<<grid_->dx<<endl;
    cout << "Polynomial  order : "<< simdata_->poly_order_  << endl;
    cout << "Runge-Kutta order : "<< simdata_->RK_order_    << endl;
    cout << "Upwind parameter  : "<< simdata_->upwind_param_<< endl;
    cout << "Poly GaussQuad order  : "<< Nquad_ << endl;
    cout << "Flux GaussQuad order  : "<< Nquad_invFlux_ << endl;
    cout << "Penalty parameter : "<< eta_penalty << endl;
    cout << "Viscous Flux scheme   : "<< simdata_->diffus_scheme_type_<<endl;
    cout <<"===============================================\n";

    return;
}

void DGSolverAdvecDiffus::InitSol(){

    register int j;
    int k=0;
    //int i;
    //double xx=0.0,qi_=0.0;

    //max_eigen_advec=0.0;  // initializing the maximum eigen value with zero

    if(simdata_->eqn_type_=="visc_burger"){  // burger's equation
        //GaussQuad quad_temp; quad_temp.setup_quadrature(8);
        //max_eigen_advec=0.0;  // initializing the maximum eigen value with zero

        for(j=0; j<grid_->Nelem; j++){
            for(k=0; k<Ndof; k++)
                Qn[j][k] = initSol_legendre_proj(j,k,quad_);

            //            for (i=0; i<quad_temp.Nq; i++){
            //                xx = 0.5 * grid_->h_j[j] * quad_temp.Gaus_pts[i] + grid_->Xc[j];
            //                qi_ = evalSolution(&Qn[j][0],xx);
            //                cout << "xx: "<<xx<<"\tqi_: "<<qi_<<endl;
            //                if(fabs(qi_)>max_eigen_advec) max_eigen_advec = fabs(qi_);
            //            }
        }
        //quad_temp.Reset_quad();

    }else if(simdata_->eqn_type_=="linear_advec_diffus"){ // linear sdvection-diffusion equation
        for(j=0; j<grid_->Nelem; j++)
            for(k=0; k<Ndof; k++)
                Qn[j][k] = initSol_legendre_proj(j,k,quad_);

        max_eigen_advec = simdata_->a_wave_;

    }else{
        FatalError_exit("Equation type is not implemented");
    }

    CalcTimeStep(); // based on maximum eigenvalues

    return;
}

double DGSolverAdvecDiffus::initSol_legendre_proj(const int &eID,
                                                  const int &basis_k,
                                                  const GaussQuad &quad_){
    int i=0,j=0,k=0;
    k=basis_k;
    j=eID;
    double xx=0.0;
    double II=0.0;
    double Qinit_=0.0;
    double Lk_=1.0;

    for (i=0; i<quad_.Nq; i++){
        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
        if(simdata_->wave_form_==3) // Decaying Burger's Turbulence, Energy model of OmerSan2016
            Qinit_ = eval_init_u_decay_burger_turb(xx);
        else
            Qinit_= eval_init_sol(xx);

        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);
        II += quad_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverAdvecDiffus::ComputeExactSolShift(){
    exact_sol_shift = (wave_speed_ * phy_time );
    return;
}

double DGSolverAdvecDiffus::ExactSol_legendre_proj(const int &eID,
                                                   const int &basis_k,
                                                   const GaussQuad &quad_){
    int i=0,j=0,k=0;
    k=basis_k;
    j=eID;
    double xx=0.0;
    double x0,x1;
    double II=0.0;
    double Qinit_=0.0;
    double Lk_=1.0;

    ComputeExactSolShift(); // updating time accurate shift

    for (i=0; i<quad_.Nq; i++){
        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i]
                + grid_->Xc[j] - exact_sol_shift;

        if(simdata_->wave_form_==0){  // single mode wave
            Qinit_ = eval_init_sol(xx);
        }else if(simdata_->wave_form_==1){  // Gaussian wave
            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);
            Qinit_= eval_init_sol(x0)+eval_init_sol(x1);
            if(x0==0 && x1==0) Qinit_ = 0.5*Qinit_;
        }else if(simdata_->wave_form_==3){
            Qinit_=0.0;
        }
        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);
        II += quad_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverAdvecDiffus::UpdateResid(double **Resid_, double **Qn_){

    register int j;
    double dQl=0.0,dQr=0.0, Ql=0.0, Qr=0.0;

    // Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    // fixme: Left and right boundary fluxes :
    j=0.0;
    Ql = evalSolution(&Qn_[grid_->Nelem-1][0], 1);
    Qr = evalSolution(&Qn_[j][0], 0);
    dQl = eval_local_du_fast(grid_->Nelem-1, &Qn_[grid_->Nelem-1][0], 1);
    dQr = eval_local_du_fast(j, &Qn_[j][0], 0);
    // Viscous Common Flux:
    int Nghost_l=1; // Needed for BR1 non-compactness
    viscflux_com[j] = Compute_common_du_flux(dQl,dQr);
    u_sol_jump[j+Nghost_l] = Compute_common_sol_jump(Ql,Qr);
    viscflux_com[grid_->Nfaces-1] = viscflux_com[j];
    u_sol_jump[grid_->Nfaces-1+Nghost_l] = u_sol_jump[j+Nghost_l];
    // Inviscid Common Flux:
    flux_com[j] = Compute_common_invflux(Ql,Qr,simdata_->a_wave_
                                         , simdata_->upwind_param_);
    flux_com[grid_->Nfaces-1] = flux_com[j];

    // Interior Faces
    for(j=1; j<grid_->Nfaces-1; j++){
        Ql = evalSolution(&Qn_[j-1][0], 1);
        Qr = evalSolution(&Qn_[j][0], 0);
        dQl = eval_local_du_fast(j-1, &Qn_[j-1][0], 1);
        dQr = eval_local_du_fast(j, &Qn_[j][0], 0);
        // Viscous Common Flux:
        viscflux_com[j] = Compute_common_du_flux(dQl,dQr);
        u_sol_jump[j+Nghost_l] = Compute_common_sol_jump(Ql,Qr);
        // Inviscid Common Flux:
        flux_com[j] = Compute_common_invflux(Ql,Qr,simdata_->a_wave_
                                             , simdata_->upwind_param_);
    }
    // Updating the additional jump values for BR1:
    u_sol_jump[0] = u_sol_jump[grid_->Nfaces-2+Nghost_l]; // for left Boundary
    u_sol_jump[grid_->Nfaces+Nghost_l] = u_sol_jump[1+Nghost_l]; // for right Boundary
    // Element loop to calculate and update the residual:
    //----------------------------------------------------
    for(j=0; j<grid_->Nelem; j++)
        UpdateResidOneCell(j, &Qn_[j][0], &Resid_[j][0]);

    return;
}

void DGSolverAdvecDiffus::UpdateResidOneCell(const int &cellid, double *q_, double *resid_){

    // General parameters:
    double inviscid_resid=0.0, viscous_resid=0.0;
    unsigned int j=cellid;
    int k=0;
    double Mkk=0.0;
    double hjj=0.0;
    double fact_=0.0;
    double Lk_p1=0.0, Lk_m1=0.0;
    double dLk_p1=0.0, dLk_m1=0.0;
    hjj = grid_->h_j[j];
    fact_ = 2.0/hjj;
    double tempC_lift  = C_lift / hjj; // to avoid division several times
    double tempC1_lift = C1_lift/hjj;

    // Preparing parameters for the Viscous Residual
    double term1=0.0,term2=0.0,term3=0.0,term4=0.0,term5=0.0;
    double du_proj_k=0.0;
    double du_flux_jp1=0.0;  // f_j+1/2
    double du_flux_jm1=0.0;  // f_j-1/2
    double u_jump_jp12 = 0.0;  // [[u]]_j+1/2
    double u_jump_jm12 = 0.0;  // [[u]]_j-1/2
    double u_jump_jp32 = 0.0;  // [[u]]_j+3/2
    double u_jump_jm32 = 0.0;  // [[u]]_j-3/2

    int Nghost_l=1; // Needed for BR1 non-compactness
    /*if(simdata_->diffus_scheme_type_=="LDG"){
        // giving that beta_e+1/2 = 1, beta_e-1/2=0
        u_jump_jm12 = 0.0;
        u_jump_jp12 = 2.0*u_sol_jump[j+1+Nghost_l];
    }else*/
    if(simdata_->diffus_scheme_type_=="SIP"
             ||simdata_->diffus_scheme_type_=="BR2"
             ||simdata_->diffus_scheme_type_=="LDG"){
        u_jump_jm12 = u_sol_jump[j+Nghost_l];
        u_jump_jp12 = u_sol_jump[j+1+Nghost_l];
    }else if(simdata_->diffus_scheme_type_=="BR1"){
        u_jump_jm12 = u_sol_jump[j+Nghost_l];
        u_jump_jp12 = u_sol_jump[j+1+Nghost_l];
        u_jump_jp32 = u_sol_jump[j+2+Nghost_l];
        u_jump_jm32 = u_sol_jump[j-1+Nghost_l];
    }

    du_flux_jm1 = viscflux_com[j];
    du_flux_jp1 = viscflux_com[j+1];

    // Preparing parameters for the InViscid Residual
    double f_proj_k=0.0;
    double flux_jp1=0.0;  // f_j+1/2
    double flux_jm1=0.0;  // f_j-1/2
    flux_jm1 = flux_com[j];
    flux_jp1 = flux_com[j+1];

    for(k=0; k<Ndof; k++){
        Mkk = 1.0/Lk_norm_squar[k];
        Lk_m1 = Lk[k][0];
        Lk_p1 = Lk[k][1];
        dLk_m1 = dLk[k][0];
        dLk_p1 = dLk[k][1];
        // Viscous Residual Computing:
        du_proj_k = eval_local_du_fluxproj_exact(j,q_,k);
        term1 = du_flux_jp1 * Lk_p1 - du_flux_jm1 *Lk_m1;
        term2 = tempC_lift*( u_jump_jp12 * Lk_p1 - u_jump_jm12 * Lk_m1); // local lift term
        term3 = tempC1_lift*( (u_jump_jp32 + u_jump_jm12 )* Lk_p1
                              -(u_jump_jp12 + u_jump_jm32 )* Lk_m1); // for BR1, 0.0 for others
        term4 = - du_proj_k ;
        if(simdata_->diffus_scheme_type_=="LDG")
            term5 = - fact_ * u_jump_jp12 * dLk_p1 ;
        else
            term5 = - 0.5 * fact_ * ( u_jump_jp12 * dLk_p1 + u_jump_jm12 * dLk_m1 );

        viscous_resid =  simdata_->thermal_diffus
                * ( term1 + term2 + term3 + term4 + term5) ;
        // Inviscid Residual Computing:
        f_proj_k = eval_local_invflux_proj_exact(q_,k);
        inviscid_resid =  - ( flux_jp1 * Lk_p1 - flux_jm1 *Lk_m1) + f_proj_k ;
        resid_[k] = fact_ * Mkk * ( viscous_resid + inviscid_resid ) ;
    }

    return;
}

double DGSolverAdvecDiffus::Compute_common_invflux(const double &Ql, const double &Qr
                                                   , const double &wave_speed
                                                   , const double &upwind_Beta_){

    double f_upw=0.0, f_cent=0.0, f_common_=0.0;
    double aa=wave_speed;
    double BB = upwind_Beta_;
    double Fl=0.0,Fr=0.0;

    if(simdata_->eqn_type_=="visc_burger"){ // burger's equation
        f_upw = Rusanov(Ql,Qr);
        Fl = 0.5 * pow(Ql,2);
        Fr = 0.5 * pow(Qr,2);

    }else if(simdata_->eqn_type_=="linear_advec_diffus"){  // linear wave equation
        f_upw = 0.5 * ( aa*(Ql+Qr) - fabs(aa) * (Qr-Ql) );
        Fl = aa * Ql;
        Fr = aa * Qr;

    }else{
        FatalError_exit("eqn type is not defined");
    }

    f_cent = 0.5 * (Fl+Fr) ;
    f_common_ = (BB * f_upw) + ((1.0-BB) * f_cent );

    return f_common_;
}

double DGSolverAdvecDiffus::Rusanov(const double &Ql, const double &Qr){

    // Now it is only working for burger's equation:
    double Lambda_max = 0.0;
    double Fl=0.0,Fr=0.0;

    Fl = 0.5 *  pow(Ql,2);
    Fr = 0.5 *  pow(Qr,2);
    Lambda_max = 0.5 * fabs(Ql+Qr);

    return  ( 0.5 * (Fl+Fr) - 0.5 * Lambda_max * (Qr-Ql) );
}

double DGSolverAdvecDiffus::eval_burgers_invflux(const double& xi_pt_, const double *q_){

    int k;
    double xx=xi_pt_;
    double Q_;

    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_basis_poly(xx,k);

    return ( 0.5 * pow(Q_,2) );
}

double DGSolverAdvecDiffus::eval_local_invflux_proj(const double *q_, const int &basis_k_){

    int k=basis_k_;
    int i;

    double II=0.0, dL_k=0.0;
    double burger_flux_=0.0;
    II=0.0;

    for (i=0; i<quad_invF_.Nq; i++){
        dL_k = eval_basis_poly_derivative(quad_invF_.Gaus_pts[i],k);
        burger_flux_ = eval_burgers_invflux(quad_invF_.Gaus_pts[i],q_);
        II += quad_invF_.Gaus_wts[i] * burger_flux_ * dL_k;
    }

    return II;
}

double DGSolverAdvecDiffus::eval_local_invflux_proj_exact(const double *q_, const int &basis_k_){

    int i;
    double II=0.0;

    if(simdata_->eqn_type_=="visc_burger"){   // burger's equation
        switch (Ndof){
        case 1:   // p0 polynomial
            switch (basis_k_) {
            case 0:
                II=0.0;
                break;
            default:
                FatalError_exit("k basis exceeds Ndof");
                break;
            }
            break;

        case 2:   // p1 polynomial
            switch (basis_k_) {
            case 0:
                II=0.0;
                break;
            case 1:
                for(i=0; i<Ndof; i++)
                    II+= pow(q_[i],2)*Lk_norm_squar[i];
                II = 0.5*II;
                break;
            default:
                FatalError_exit("k basis exceeds Ndof");
                break;
            }
            break;

        case 3:   // p2 polynomial
            switch (basis_k_) {
            case 0:
                II=0.0;
                break;
            case 1:
                for(i=0; i<Ndof; i++)
                    II+= pow(q_[i],2)*Lk_norm_squar[i];
                II = 0.5*II;
                break;
            case 2:
                for(i=0; i<Ndof-1; i++){
                    II+= q_[i]*q_[i+1]*Lk_norm_squar[i]*Lk_norm_squar[i+1]*(i+1);
                }
                II = (3./2.) * II;
                break;
            default:
                FatalError_exit("k basis exceeds Ndof");
                break;
            }
            break;

        case 4:      // p3 polynomial
            switch (basis_k_) {
            case 0:
                II=0.0;
                break;
            case 1:
                for(i=0; i<Ndof; i++)
                    II+= pow(q_[i],2)*Lk_norm_squar[i];
                II = 0.5*II;
                break;
            case 2:
                for(i=0; i<Ndof-1; i++){
                    II+= q_[i]*q_[i+1]*Lk_norm_squar[i]*Lk_norm_squar[i+1]*(i+1);
                }
                II = (3./2.) * II;
                break;
            case 3:
                II = pow(q_[0],2)+pow(q_[1],2)
                        +(17./35.)*pow(q_[2],2)
                        +(pow(q_[3],2)/3.)
                        +2.*q_[0]*q_[2]
                        +(6./7.)*q_[1]*q_[3];
                break;
            default:
                FatalError_exit("k basis exceeds Ndof");
                break;
            }
            break;

        case 5:  // p4 polynomial
            switch (basis_k_) {
            case 0:
                II=0.0;
                break;
            case 1:
                for(i=0; i<Ndof; i++)
                    II+= pow(q_[i],2)*Lk_norm_squar[i];
                II = 0.5*II;
                break;
            case 2:
                for(i=0; i<Ndof-1; i++){
                    II+= q_[i]*q_[i+1]*Lk_norm_squar[i]*Lk_norm_squar[i+1]*(i+1);
                }
                II = (3./2.) * II;
                break;
            case 3:
                II = pow(q_[0],2)+pow(q_[1],2)
                        +(17./35.)*pow(q_[2],2)
                        +(pow(q_[3],2)/3.)
                        +(59./231.)*pow(q_[4],2)
                        +2.*q_[0]*q_[2]
                        +(6./7.)*q_[1]*q_[3]
                        +(4./7.)*q_[2]*q_[4];
                break;
            case 4:
                II = 2.*(q_[0]*q_[1] + q_[1]*q_[2] + q_[0]*q_[3])
                        +(22./21.)*q_[2]*q_[3]
                        +(8./9.)*q_[1]*q_[4]
                        +(172./231.)*q_[3]*q_[4];
                break;
            default:
                FatalError_exit("k basis exceeds Ndof");
                break;
            }
            break;

        case 6: // p5 polynomial
            switch (basis_k_) {
            case 0:
                II=0.0;
                break;
            case 1:
                for(i=0; i<Ndof; i++)
                    II+= pow(q_[i],2)*Lk_norm_squar[i];
                II = 0.5*II;
                break;
            case 2:
                for(i=0; i<Ndof-1; i++){
                    II+= q_[i]*q_[i+1]*Lk_norm_squar[i]*Lk_norm_squar[i+1]*(i+1);
                }
                II = (3./2.) * II;
                break;
            case 3:
                II = pow(q_[0],2)+pow(q_[1],2)
                        +(17./35.)*pow(q_[2],2)
                        +(pow(q_[3],2)/3.)
                        +(59./231.)*pow(q_[4],2)
                        +(89./429.) * pow(q_[5],2)
                        +2.*q_[0]*q_[2]
                        +(6./7.)*q_[1]*q_[3]
                        +(4./7.)*q_[2]*q_[4]
                        +(100./231.)*q_[3]*q_[5];
                break;
            case 4:
                II = 2.*(q_[0]*q_[1] + q_[1]*q_[2] + q_[0]*q_[3])
                        +(22./21.)*q_[2]*q_[3]
                        +(8./9.)*q_[1]*q_[4]
                        +(172./231.)*q_[3]*q_[4]
                        +(20./33.)*q_[2]*q_[5]
                        +(250./429.)*q_[4]*q_[5];
                break;
            case 5:
                II = pow(q_[0],2)+pow(q_[1],2)+pow(q_[2],2)
                        +(131./231.)*pow(q_[3],2)
                        +(179./429.)*pow(q_[4],2)
                        +(1./3.)*pow(q_[5],2)
                        +2.*(q_[0]*q_[2]+q_[1]*q_[3]+q_[0]*q_[4])
                        +(12./11.)*q_[2]*q_[4]
                        +(10./11.)*q_[1]*q_[5]
                        +(340./429.)*q_[3]*q_[5];
                break;
            default:
                FatalError_exit("k basis exceeds Ndof");
                break;
            }
            break;

        default:
            FatalError_exit("Polynomial order is not implemented here");
            break;
        }

    }else if(simdata_->eqn_type_=="linear_advec_diffus"){ //linear advection-diffusion equation
        switch (basis_k_) {
        case 0:
            II= 0.0;
            break;
        case 1:
            II= ( 2.0 * simdata_->a_wave_ * q_[0] );
            break;
        case 2:
            II= ( 2.0 * simdata_->a_wave_ * q_[1] );
            break;
        case 3:
            II= ( 2.0 * simdata_->a_wave_ * ( q_[0] + q_[2] ) );
            break;
        case 4:
            II= ( 2.0 * simdata_->a_wave_ * ( q_[1] + q_[3] ) );
            break;
        case 5:
            II= ( 2.0 * simdata_->a_wave_ * ( q_[0] + q_[2] + q_[4] ) );
            break;
        default:
            FatalError_exit("Polynomial order is not implemented here");
            break;
        }
    }else{
        FatalError_exit("eqn type is not implemented");
    }

    return II;
}

double DGSolverAdvecDiffus::eval_local_du(const int eID,
                                          const double *q_
                                          , const double& xi_pt){
    double fact_ = 2.0/grid_->h_j[eID];
    double II=0.0;
    int k;

    for(k=0; k<Ndof; k++)
        II += q_[k] * eval_basis_poly_derivative(xi_pt, k);

    return (fact_*II);
}

double DGSolverAdvecDiffus::eval_local_du_fast(const int eID,
                                               const double *q_
                                               , const int& i_pos_){
    double fact_ = 2.0/grid_->h_j[eID];
    double II=0.0;
    int k;

    for(k=0; k<Ndof; k++)
        II += q_[k] * dLk[k][i_pos_];

    return (fact_*II);
}

double DGSolverAdvecDiffus::Compute_common_sol_jump(const double &ul_
                                                    , const double &ur_){
    return (ur_ - ul_);
}

double DGSolverAdvecDiffus::Compute_common_du_flux(const double& dul_
                                                   , const double& dur_){
    if(simdata_->diffus_scheme_type_=="SIP"
            || simdata_->diffus_scheme_type_=="BR1"
            || simdata_->diffus_scheme_type_=="BR2"){
        return ( 0.5 * ( dul_ + dur_) );
    }else if(simdata_->diffus_scheme_type_=="LDG"){
        return dul_;
    }else{
        _notImplemented("diffusion Scheme type");
        return 0.0;
    }
}

void DGSolverAdvecDiffus::Compute_vertex_sol(){

    register int i;
    double Ql=0.0, Qr=0.0;

    i=0;
    Ql = evalSolution(&Qn[grid_->Nelem-1][0],1.0);
    Qr = evalSolution(&Qn[i][0],-1.0);
    Qv[i] = 0.5 * ( Ql + Qr);
    Qv[grid_->Nfaces-1] = Qv[i];

    for(i=1; i<grid_->Nfaces-1; i++){

        Ql = evalSolution(&Qn[i-1][0],1.0);
        Qr = evalSolution(&Qn[i][0],1.0);
        Qv[i] = 0.5 * ( Ql + Qr);
    }

    return;
}

void DGSolverAdvecDiffus::Compute_exact_vertex_sol(){
    register int j;
    double xx=0.0;
    double x0,x1;

    ComputeExactSolShift(); // updating the shift for the currnet time

    for(j=0; j<grid_->N_exact_ppts; j++){
        xx = grid_->x_exact_ppts[j]- exact_sol_shift;

        if(simdata_->wave_form_==0){  // single wave mode
            Q_exact[j] = eval_init_sol(xx);
        }else if(simdata_->wave_form_==1){  // Gaussian wave
            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);
            if(x0==0 && x1==0)
                Q_exact[j] = 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
            else
                Q_exact[j] = (eval_init_sol(x0)+ eval_init_sol(x1));

        }else if(simdata_->wave_form_==3){ // burger's decay turb
            Q_exact[j] = eval_init_u_decay_burger_turb(xx);

        }else{
            FatalError_exit("Wave form is not implemented");
        }
    }

    return;
}

void DGSolverAdvecDiffus::Compute_projected_exact_sol(){
    register int j; int k=0;
    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qex_proj[j][k] = ExactSol_legendre_proj(j,k,quad_);
    return;
}

double DGSolverAdvecDiffus::eval_init_sol(const double& xx){

    if(simdata_->wave_form_==0){  // u(x,0) = A * sin ( (f * PI * x) / L + phy ) + C
        double argum_ = (simdata_->wave_freq_*PI*xx)
                / fabs(wave_length_) + simdata_->wave_shift ;
        double wave_value_ = simdata_->wave_amp_ * sin(argum_)
                + simdata_->wave_const;
        return wave_value_;
    }else if(simdata_->wave_form_==1){
        double wave_value_ = simdata_->Gaussian_amp_
                * exp(-simdata_->Gaussian_exponent_*pow(xx,2)) ;
        return wave_value_;
    }else{
        _notImplemented("Wave form is not implemented");
    }
}

double DGSolverAdvecDiffus::eval_init_u_decay_burger_turb(const double& xx_){

    register int i;
    double u_=0.;
    double dk_,k_max_, E_, k_, epsi_;
    k_max_ = simdata_->max_wave_no_;
    dk_ = 1.0;
    int n_pts_=k_max_/dk_;

    for(i=0; i<n_pts_; i++){
        k_ = simdata_->k_wave_no_[i];
        epsi_= simdata_->epsi_phase_[i];
        E_= simdata_->energy_spect_[i];
        u_ += sqrt(2.*E_ ) * cos (k_ * xx_ + 2.*PI*epsi_);
    }

    return (u_ + simdata_->velocity_mean_ );
}

double DGSolverAdvecDiffus::eval_exact_sol(double &xx){

    ComputeExactSolShift(); // updating the shift
    xx = xx - exact_sol_shift;
    if(simdata_->wave_form_==0){
        return eval_init_sol(xx);
    }else if(simdata_->wave_form_==1){
        double x0 = xx - wave_length_*floor(xx/wave_length_);
        double x1 = xx + wave_length_*floor(xx/-wave_length_);
        if(x0==0 && x1==0)
            return 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
        else
            return (eval_init_sol(x0)+ eval_init_sol(x1));
    }else if(simdata_->wave_form_==3){
        return 0.0;
    }else{
        _notImplemented("Wave form is not implemented");
        return 0.0;
    }
}

double DGSolverAdvecDiffus::eval_basis_poly(const double& xi_, const int& basis_k_){

    double Lk_=0.0;

    switch (basis_k_) {
    case 0:
        Lk_ = 1.0;  break;
    case 1:
        Lk_ = xi_;  break;
    case 2:
        Lk_ = 0.5 * (3.0 * xi_*xi_ -1.0);  break;
    case 3:
        Lk_ = 0.5 * (5.0 * pow(xi_,3) - 3.0 * xi_);  break;
    case 4:
        Lk_ = (35. * pow(xi_,4) - 30. * xi_*xi_ + 3. ) /8. ;  break;
    case 5:
        Lk_ = (63.0 * pow(xi_,5) - 70.0 * pow(xi_,3) + 15.0 * xi_ ) /8. ;  break;
    default:
        FatalError_exit("k basis exceeds poly order");
        break;
    }

    return Lk_;
}

double DGSolverAdvecDiffus::eval_basis_poly_derivative(const double& xi_
                                                       , const int& basis_k_){
    double dLk_=1000.0;

    switch (basis_k_) {
    case 0: dLk_= 0.0;  break;
    case 1: dLk_= 1.0;  break;
    case 2: dLk_= 3.0 * xi_;  break;
    case 3: dLk_= 0.5 * (15.0 * pow(xi_,2) - 3.0);  break;
    case 4: dLk_= 0.5 * (35.0 * pow(xi_,3) - 15.0 * xi_ );  break;
    case 5: dLk_= ( 315.0 * pow(xi_,4) - 210.0 * pow(xi_,2) + 15.0 ) / 8. ;  break;
    default:
        FatalError_exit("k basis exceeds poly order");
        break;
    }

    return dLk_;
}

double DGSolverAdvecDiffus::eval_basis_norm_squared(const int &basis_k_){
    // this is norm^2
    return 2./(2*basis_k_+1);
}

double DGSolverAdvecDiffus::evalSolution(const double* q_, const double& xi_pt_){

    double xx=xi_pt_;
    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_basis_poly(xx,k);

    return Q_;
}

double DGSolverAdvecDiffus::evalSolution(const double* q_, const int& i_pos_){

    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * Lk[k][i_pos_];   // i_pos_ = 0 if xi =-1, i_pos_ =1 if xi =1

    return Q_;
}

double DGSolverAdvecDiffus::eval_local_du_fluxproj(const int eID, const double *q_
                                                   , const int &basis_k_){
    int k=basis_k_;
    int i;
    double II=0.0, dL_k=0.0;
    double du_flux_=0.0;

    for (i=0; i<quad_viscF_.Nq; i++){
        dL_k = eval_basis_poly_derivative(quad_viscF_.Gaus_pts[i],k);
        du_flux_ = eval_local_du(eID, q_, quad_viscF_.Gaus_pts[i]);
        II += quad_viscF_.Gaus_wts[i] * du_flux_ * dL_k;
    }

    return II;
}

double DGSolverAdvecDiffus::eval_local_du_fluxproj_exact(const int eID, const double *q_
                                                         , const int &basis_k_){
    int i;
    double II=0.0;
    double fact_ = 2.0/grid_->h_j[eID];

    switch(basis_k_){
    case 0: // p0 -> p5
        II = 0.0;
        break;
    case 1:  // p1 -> p5
        for(i=1; i<Ndof; i+=2)
            II+= q_[i];
        II = 2.0* fact_ * II;
        break;
    case 2:  // p2 -> p5
        for(i=2; i<Ndof; i+=2)
            II+= q_[i];
        II = 6.0 * fact_ * II;
        break;
    case 3:  // p3 -> p5
        switch(Ndof){
        case 4: // p3
        case 5: // p4
            II = 2.0 *fact_ *(q_[1]+6.0*q_[3]);
            break;
        case 6: // p5
            II = 2.0 *fact_ *(q_[1]+6.0*(q_[3]+q_[5]));
            break;
        default:
            FatalError_exit("Ndof is not correct");
            break;
        }
        break;
    case 4: // p4, p5
        II = fact_ * (6.0 * q_[2] + 20.0 * q_[4]);
        break;
    case 5: // p5
        II = 2.0 *fact_ *(q_[1]+6.0*q_[3]+15.0*q_[5]);
        break;
    default:
        FatalError_exit("k basis exceeds poly order");
        break;
    }

    return II;
}

double DGSolverAdvecDiffus::ComputePolyError(){

    register int j; int i;
    double xx=0.0,L2_error=0.0,elem_error=0.0,II=0.0,q_ex,q_n;

    for(j=0; j<grid_->Nelem; j++){
        elem_error=0.0;
        for(i=0; i<quad_.Nq; i++) {
            xx= 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i]
                    + grid_->Xc[j] - exact_sol_shift;
            q_ex = eval_init_sol(xx);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            elem_error += quad_.Gaus_wts[i] * pow((q_ex - q_n),2);
        }
        II += (0.5 * grid_->h_j[j] * elem_error) ;
    }
    L2_error = sqrt(II/(grid_->xf-grid_->x0));

    return L2_error;
}

double DGSolverAdvecDiffus::L1_error_nodal_gausspts_proj(){

    register int j; int i;
    double L1_error=0.0,II=0.0,q_ex,q_n;
    II=0.0;

    for(j=0; j<grid_->Nelem; j++){
        for(i=0; i<quad_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j][0], quad_.Gaus_pts[i]);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            II += fabs(q_ex - q_n);
        }
    }
    L1_error = II/(Ndof*grid_->Nelem);

    return L1_error;
}

double DGSolverAdvecDiffus::L2_error_nodal_gausspts_proj(){

    register int j; int i;

    double L2_error=0.0,II=0.0,q_ex,q_n;

    II=0.0;

    for(j=0; j<grid_->Nelem; j++){
        for(i=0; i<quad_.Nq; i++) {

            q_ex = evalSolution(&Qex_proj[j][0], quad_.Gaus_pts[i]);

            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);

            II += pow((q_ex - q_n),2);
        }
    }

    L2_error = sqrt(II/(Ndof*grid_->Nelem));

    return L2_error;
}

double DGSolverAdvecDiffus::L1_error_nodal_gausspts(){

    register int j; int i;

    double L1_error=0.0,II=0.0,q_ex,q_n,xx=0.0;

    II=0.0;

    for(j=0; j<grid_->Nelem; j++){
        for(i=0; i<quad_.Nq; i++) {
            xx =  0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
            q_ex = eval_exact_sol(xx);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            II += fabs(q_ex - q_n);
        }
    }

    L1_error = II/(Ndof*grid_->Nelem);

    return L1_error;
}

double DGSolverAdvecDiffus::L2_error_nodal_gausspts(){

    register int j; int i;

    double L2_error=0.0,II=0.0,q_ex,q_n,xx=0.0;

    II=0.0;

    for(j=0; j<grid_->Nelem; j++){
        for(i=0; i<quad_.Nq; i++) {
            xx =  0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
            q_ex = eval_exact_sol(xx);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            II += pow((q_ex - q_n),2);
        }
    }

    L2_error = sqrt(II/(Ndof*grid_->Nelem));
    return L2_error;
}

double DGSolverAdvecDiffus::L1_error_projected_sol(){

    register int j; int i;

    double L1_error=0.0,elem_error=0.0,II=0.0,q_ex,q_n;

    for(j=0; j<grid_->Nelem; j++){

        elem_error=0.0;
        for(i=0; i<quad_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j][0], quad_.Gaus_pts[i]);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            elem_error += quad_.Gaus_wts[i] * fabs(q_ex - q_n);
        }

        II += (0.5 * grid_->h_j[j] * elem_error) ;
    }

    L1_error = II/(grid_->xf-grid_->x0);

    return L1_error;
}

double DGSolverAdvecDiffus::L2_error_projected_sol(){

    register int j; int i;
    double L2_error=0.0,elem_error=0.0,II=0.0,q_ex,q_n;

    for(j=0; j<grid_->Nelem; j++){

        elem_error=0.0;
        for(i=0; i<quad_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j][0], quad_.Gaus_pts[i]);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            elem_error += quad_.Gaus_wts[i] * pow((q_ex - q_n),2);
        }
        II += (0.5 * grid_->h_j[j] * elem_error) ;
    }
    L2_error = sqrt(II/(grid_->xf-grid_->x0));

    return L2_error;
}

double DGSolverAdvecDiffus::L1_error_average_sol(){

    register int j;

    double L1_error=0.0,error=0.0;

    for(j=0; j<grid_->Nelem; j++){

        error += fabs(Qex_proj[j][0] - Qn[j][0]);
    }

    L1_error = error/grid_->Nelem;

    return L1_error;
}

double DGSolverAdvecDiffus::L2_error_average_sol(){

    register int j;

    double L2_error=0.0,error=0.0;

    for(j=0; j<grid_->Nelem; j++){

        error += pow((Qex_proj[j][0] - Qn[j][0]),2);
    }

    L2_error = sqrt(error/grid_->Nelem);

    return L2_error;
}

void DGSolverAdvecDiffus::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->wave_form_!=3){
        if(simdata_->Sim_mode=="error_analysis_dt"){

            sprintf(fname,"%snodal/u_cont_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                    ,simdata_->case_postproc_dir
                    ,grid_->Nelem
                    ,time_step
                    ,eta_penalty
                    ,simdata_->Nperiods);

            FILE* sol_out=fopen(fname,"w");

            for(j=0; j<grid_->Nfaces; j++)
                fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qv[j]);

            fclose(sol_out);
            emptyarray(fname);

        } else{

            sprintf(fname,"%snodal/u_cont_N%d_CFL%1.3e_Beta%1.2f_Eps%1.2f_%1.3fT.dat"
                    ,simdata_->case_postproc_dir
                    ,grid_->Nelem
                    ,CFL
                    ,simdata_->upwind_param_
                    ,eta_penalty
                    ,simdata_->Nperiods);

            FILE* sol_out=fopen(fname,"w");

            for(j=0; j<grid_->Nfaces; j++)
                fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qv[j]);

            fclose(sol_out);
            emptyarray(fname);
        }
    }

    fname = new char[100];

    sprintf(fname,"%snodal/u_cont_exact_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,simdata_->Nperiods);

    FILE* sol_out1=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out1, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact[j]);

    fclose(sol_out1);

    emptyarray(fname);

    return;
}

void DGSolverAdvecDiffus::print_average_sol(){

    register int j;

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%saver/u_aver_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++)
            fprintf(sol_out, "%2.10e %2.10e %2.10e\n"
                    ,grid_->Xc[j], Qex_proj[j][0], Qn[j][0]);

        fclose(sol_out);
        emptyarray(fname);

    }else{

        sprintf(fname,"%saver/u_aver_N%d_CFL%1.3e_Beta%1.2f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++)
            fprintf(sol_out, "%2.10e %2.10e %2.10e\n"
                    ,grid_->Xc[j], Qex_proj[j][0], Qn[j][0]);

        fclose(sol_out);
        emptyarray(fname);
    }


    return;
}

void DGSolverAdvecDiffus::dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                                      ,double& L1_aver_sol_,double& L2_aver_sol_
                                      ,double& L1_nodal_gausspts, double& L2_nodal_gausspts){

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="error_analysis_CFL"){
        sprintf(fname,"%serrors/errors_CFL%1.3f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,CFL
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%d %2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,grid_->Nelem
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);
        emptyarray(fname);

    }else if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%serrors/errors_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"w");

        fprintf(solerror_out, "%2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);
        emptyarray(fname);

        // Dumping all errors in one file as a function of dt:
        //--------------------------------------------------------
        fname = new char[100];

        sprintf(fname,"%serrors/errors_N%d_alldt_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,eta_penalty
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%1.7e %2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,time_step
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);

        emptyarray(fname);

        // Dumping all errors in one file as a function of Nelem:
        //--------------------------------------------------------
        fname = new char[100];

        sprintf(fname,"%serrors/errors_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,time_step
                ,eta_penalty
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%d %2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,grid_->Nelem
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);

        emptyarray(fname);

    }else if( simdata_->Sim_mode=="test" || simdata_->Sim_mode=="normal" ){
        sprintf(fname,"%serrors/errors_N%d_CFL%1.3e_Beta%1.2f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"w");

        fprintf(solerror_out, "%2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);

        emptyarray(fname);

    }else if(simdata_->Sim_mode=="error_analysis_Eps"){

        sprintf(fname,"%serrors/errors_N%d_CFL%1.3f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"w");

        fprintf(solerror_out, "%2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);

        emptyarray(fname);

        // Dumping all errors in one file as a function of beta:
        //--------------------------------------------------------
        fname = new char[100];

        sprintf(fname,"%serrors/errors_N%d_CFL%1.3f_allEps_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%1.2f %2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,simdata_->upwind_param_
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);
        emptyarray(fname);
    }

    return;
}

void DGSolverAdvecDiffus::dump_discont_sol(){

    register int j; int k;

    double xx=0.0,qq=0.0;

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%snodal/u_disc_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++){

            for(k=0; k<grid_->N_xi_disc_ppts; k++) {

                qq = evalSolution(&Qn[j][0],grid_->xi_disc[k]);

                xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                        + grid_->Xc[j];

                fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
            }
        }

        fclose(sol_out);
        emptyarray(fname);

    }else {

        sprintf(fname,"%snodal/u_disc_N%d_CFL%1.3e_Beta%1.2f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,eta_penalty
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++){

            for(k=0; k<grid_->N_xi_disc_ppts; k++) {

                qq = evalSolution(&Qn[j][0],grid_->xi_disc[k]);

                xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                        + grid_->Xc[j];

                fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
            }
        }

        fclose(sol_out);
        emptyarray(fname);
    }

    fname = new char[100];

    sprintf(fname,"%snodal/u_disc_exact_N%d_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,simdata_->Nperiods);

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nelem; j++){

        for(k=0; k<grid_->N_xi_disc_ppts; k++) {

            qq = evalSolution(&Qex_proj[j][0],grid_->xi_disc[k]);

            xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                    + grid_->Xc[j];

            fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
        }
    }

    fclose(sol_out);
    emptyarray(fname);

    return;
}

void DGSolverAdvecDiffus::dump_timeaccurate_sol(){

    register int j; int k;

    double xx=0.0,qq=0.0;

    char *fname=nullptr;
    fname = new char[250];

    // Dump continuous data on uniform points:
    //-------------------------------------------
    compute_uniform_cont_sol(); // compute time accurate solution

    if(simdata_->Sim_mode=="CFL_const"
            || simdata_->Sim_mode =="error_analysis_CFL"
            || simdata_->Sim_mode =="test"
            || simdata_->Sim_mode =="normal"){
        // Dump time accurate continuous equally spaced solution data:
        sprintf(fname,"%stime_data/u_cont_N%d_CFL%1.4f_Beta%1.2f_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        // Dump time accurate continuous equally spaced solution data:
        sprintf(fname,"%stime_data/u_cont_N%d_dt%1.3e_Beta%1.2f_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
                ,eta_penalty
                ,phy_time);

    }

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->N_uniform_pts; j++)
        fprintf(sol_out,"%2.10e %2.10e\n"
                ,grid_->x_unifrom_pts[j],Q_cont_sol[j]);

    fclose(sol_out);
    emptyarray(fname);

    // Dump time accurate Discontinuous data:
    fname = new char[250];
    if(simdata_->Sim_mode=="CFL_const"
            || simdata_->Sim_mode =="error_analysis_CFL"
            || simdata_->Sim_mode =="test"
            || simdata_->Sim_mode =="normal"){
        sprintf(fname,"%stime_data/u_disc_N%d_CFL%1.4f_Beta%1.2f_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        sprintf(fname,"%stime_data/u_disc_N%d_dt%1.3e_Beta%1.2f_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
                ,eta_penalty
                ,phy_time);
    }

    sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<grid_->N_xi_disc_ppts; k++) {
            qq = evalSolution(&Qn[j][0],grid_->xi_disc[k]);
            xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                    + grid_->Xc[j];
            fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
        }

    fclose(sol_out);
    emptyarray(fname);

    return;
}

void DGSolverAdvecDiffus::compute_uniform_cont_sol(){

    // Dump continuous data on uniform points:
    //-------------------------------------------
    register int j; int k;

    int count_=0; //continuous points counter

    // For element zero:
    k = simdata_->N_uniform_pts_per_elem_-1;
    j= grid_->Nelem-1;
    Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]); // a potential bug, fixed by changing j-1 to j

    k=0; j=0;
    Q_cont_sol[count_] += evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
    Q_cont_sol[count_] = 0.5*Q_cont_sol[count_];
    count_++;

    for(k=1; k<simdata_->N_uniform_pts_per_elem_-1; k++) {
        Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
        count_++;
    }

    for(j=1; j<grid_->Nelem; j++){

        k=0;
        Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
        k = simdata_->N_uniform_pts_per_elem_-1;
        Q_cont_sol[count_] += evalSolution(&Qn[j-1][0],grid_->xi_uniform[k]);

        Q_cont_sol[count_] = 0.5*Q_cont_sol[count_];
        count_++;

        for(k=1; k<simdata_->N_uniform_pts_per_elem_-1; k++) {
            Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
            count_++;
        }
    }
    Q_cont_sol[count_] = Q_cont_sol[0];

    //printf("\n Count: %d \t N_uniform: %d\n", count_, grid_->N_uniform_pts);

    return;
}

double DGSolverAdvecDiffus::compute_totalVariation(){

    register int i;

    double TV_=0.0;
    max_eigen_advec = 0.0;

    for(i=1; i<grid_->N_uniform_pts; i++){
        TV_ += fabs(Q_cont_sol[i]-Q_cont_sol[i-1]);
        if(fabs(Q_cont_sol[i-1])>max_eigen_advec)
            max_eigen_advec=fabs(Q_cont_sol[i-1]);
    }

    return TV_;
}























