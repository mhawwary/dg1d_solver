#include "DGSolverAdvecDiffus.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

DGSolverAdvecDiffus::~DGSolverAdvecDiffus(void){

    Reset_solver();
}

void DGSolverAdvecDiffus::setup_solver(GridData& meshdata_, SimData& osimdata_){

    grid_ = &meshdata_;
    simdata_ = &osimdata_;

    Ndof = simdata_->poly_order_+1;

    Qn    =  new double* [grid_->Nelem];

    Qex_proj = new double*[grid_->Nelem];

    register int i;

    for(i=0; i<grid_->Nelem; i++){
        Qn[i]       = new double[Ndof];
        Qex_proj[i] = new double[Ndof];
    }


    Q_exact = new double[grid_->N_exact_ppts];
    Qv = new double[grid_->Nfaces];

    flux_com = new double[grid_->Nfaces];
    viscflux_com = new double[grid_->Nfaces];
    u_sol_jump = new double[grid_->Nfaces];

    if(simdata_->diffus_scheme_type_=="SIP")
        r_lift = pow(simdata_->poly_order_,2);
    else if(simdata_->diffus_scheme_type_=="BR2")
        r_lift = pow((simdata_->poly_order_+1),2)/2.0;

    e_penalty = simdata_->penalty_param_;

    SetPhyTime(simdata_->t_init_);

    CalcTimeStep();

    ComputeExactSolShift();
    Compute_projected_exact_sol();
    Compute_exact_vertex_sol();

    // Screen Output of input and simulation parameters:
    cout <<"\n===============================================\n";
    cout << "CFL no.        : "<<CFL<<endl;
    cout << "time step, dt  : "<<time_step<<endl;
    cout << "last_time_step: "<<last_time_step<<endl;
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
    cout << "Penalty parameter : "<< e_penalty << endl;
    cout <<"===============================================\n";

    return;
}

void DGSolverAdvecDiffus::Reset_solver(){

    emptyarray(xi);

    emptyarray(grid_->Nelem,Qn);
    emptyarray(Q_exact);
    emptyarray(flux_com);
    emptyarray(viscflux_com);
    emptyarray(u_sol_jump);
    emptyarray(Qv);
    emptyarray(grid_->Nelem,Qex_proj);

    grid_->Reset_();


    return;
}


// Solver functions
//-------------------------------------------

void DGSolverAdvecDiffus::CalcTimeStep(){

    double dx = grid_->dx;
    double dx2 = pow(dx,2);
    double radius_advec_=0.0, radius_diffus_ =0.0;

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;
    radius_advec_ =  simdata_->a_wave_ / dx  ;
    radius_diffus_ = simdata_->thermal_diffus / dx2 ;

    if(simdata_->calc_dt_flag==1){

        CFL = simdata_->CFL_;

        time_step = CFL / ( radius_advec_ + radius_diffus_ );
        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){

        time_step = simdata_->dt_;
        last_time_step = time_step;

        CFL = time_step * ( radius_advec_ + radius_diffus_ );
        simdata_->CFL_ = CFL;

    }else {

        FatalError_exit("Wrong Calc_dt_flag");
    }

    // Determining end of simulation parameters:
    //----------------------------------------------------

    if(simdata_->end_of_sim_flag_==0){

        simdata_->t_end_ = simdata_->Nperiods * T_period;

        simdata_->maxIter_ = (int) ceil(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step) > simdata_->t_end_ ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step) < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==1){

        simdata_->Nperiods = simdata_->t_end_/T_period;
        simdata_->maxIter_ = (int) ceil(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step) > simdata_->t_end_ ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step) < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==2){

        simdata_->t_end_ = simdata_->maxIter_ * time_step;
        simdata_->Nperiods = simdata_->t_end_/T_period;

    }else{
        FatalError_exit("Wrong end_of_simulation_flag");
    }

    return;
}

void DGSolverAdvecDiffus::InitSol(){

    register int j;

    int k=0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    for(j=0; j<grid_->Nelem; j++){

        for(k=0; k<Ndof; k++){

            Qn[j][k] = initSol_legendre_proj(j,k,quad_);
        }
    }

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
    double Lk_norm=1.0;  // norm^2

    II=0.0;

    for (i=0; i<quad_.Nq; i++){

        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
        Qinit_= eval_init_sol(xx);

        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);

        II += quad_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }

    Lk_norm = eval_basis_norm_squared(k);
    II = II / Lk_norm ;

    return II;
}

void DGSolverAdvecDiffus::ComputeExactSolShift(){

    // Preparing shift information:
    //-------------------------------
    double a=0.;

    wave_length_ = grid_->xf - grid_->x0 ;
    a = simdata_->a_wave_;
    exact_sol_shift = (a * simdata_->t_end_ );

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
    double Lk_norm=1.0;  // norm^2

    II=0.0;

    for (i=0; i<quad_.Nq; i++){

        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j] - exact_sol_shift;

        if(simdata_->wave_form_==0){

            Qinit_ = eval_init_sol(xx);

        }else if(simdata_->wave_form_==1){  // Gaussian wave

            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);

            Qinit_= eval_init_sol(x0)+eval_init_sol(x1);

            if(x0==0 && x1==0) Qinit_ = 0.5*Qinit_;
        }

        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);

        II += quad_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }

    Lk_norm = eval_basis_norm_squared(k);
    II = II / Lk_norm ;

    return II;
}

void DGSolverAdvecDiffus::UpdateResid(double **Resid_, double **Qn_){

    register int j;

    // Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    double dQl=0.0,dQr=0.0, Ql=0.0, Qr=0.0;

    // fixme: Left and right boundary fluxes :

    j=0.0;

    Ql = evalSolution(&Qn_[grid_->Nelem-1][0], 1.0);
    Qr = evalSolution(&Qn_[j][0], -1.0);
    dQl = eval_local_du(grid_->Nelem-1, &Qn_[grid_->Nelem-1][0], 1.0);
    dQr = eval_local_du(j, &Qn_[j][0], -1.0);

    // Viscous Common Flux:
    viscflux_com[j] = Compute_common_du_flux(dQl,dQr);
    u_sol_jump[j] = Compute_common_sol_jump(Ql,Qr);

    viscflux_com[grid_->Nfaces-1] = viscflux_com[j];
    u_sol_jump[grid_->Nfaces-1] = u_sol_jump[j];

    // Inviscid Common Flux:
    flux_com[j] = Compute_common_flux(Ql,Qr,simdata_->a_wave_
                                      , simdata_->upwind_param_);
    flux_com[grid_->Nfaces-1] = flux_com[j];

    for(j=1; j<grid_->Nfaces-1; j++){

        Ql = evalSolution(&Qn_[j-1][0], 1.0);
        Qr = evalSolution(&Qn_[j][0], -1.0);
        dQl = eval_local_du(j-1, &Qn_[j-1][0], 1.0);
        dQr = eval_local_du(j, &Qn_[j][0], -1.0);

        // Viscous Common Flux:
        viscflux_com[j] = Compute_common_du_flux(dQl,dQr);
        u_sol_jump[j] = Compute_common_sol_jump(Ql,Qr);
        // Inviscid Common Flux:
        flux_com[j] = Compute_common_flux(Ql,Qr,simdata_->a_wave_
                                          , simdata_->upwind_param_);
    }


    // Element loop to calculate and update the residual:
    //----------------------------------------------------

    for(j=0; j<grid_->Nelem; j++){
         UpdateResidOneCell(j, &Qn_[j][0], &Resid_[j][0]);
    }

    return;
}

void DGSolverAdvecDiffus::UpdateResidOneCell(const int &cellid, double *q_, double *resid_){

    // General parameters:
    double inviscid_resid=0.0, viscous_resid=0.0;

    unsigned int j=cellid;

    int k=0;

    double Mkk=0.0,mkk=0.0;
    double hjj=0.0;
    double fact_=0.0;
    double Lk_p1=0.0, Lk_m1=0.0;
    double dLk_p1=0.0, dLk_m1=0.0;

    hjj = grid_->h_j[j];
    fact_ = 2.0/hjj;

    // Preparing parameters for the Viscous Residual
    double term1=0.0,term2=0.0,term3=0.0,term4=0.0;

    double du_proj_k=0.0;
    double du_flux_jp1=0.0;  // f_j+1/2
    double du_flux_jm1=0.0;  // f_j-1/2
    double u_sol_jp1 = 0.0;
    double u_sol_jm1 = 0.0;

    u_sol_jm1 = u_sol_jump[j];
    u_sol_jp1 = u_sol_jump[j+1];
    du_flux_jm1 = viscflux_com[j];
    du_flux_jp1 = viscflux_com[j+1];

    // Preparing parameters for the InViscid Residual
    double f_proj_k=0.0;
    double flux_jp1=0.0;  // f_j+1/2
    double flux_jm1=0.0;  // f_j-1/2

    flux_jm1 = flux_com[j];
    flux_jp1 = flux_com[j+1];

    for(k=0; k<Ndof; k++){

        mkk = eval_basis_norm_squared(k);
        Mkk = 1.0/mkk;
        Lk_m1 = eval_basis_poly(-1.0,k);
        Lk_p1 = eval_basis_poly( 1.0,k);
        dLk_m1 = eval_basis_poly_derivative(-1.0,k);
        dLk_p1 = eval_basis_poly_derivative(1.0,k);

        // Viscous Residual Computing:
        du_proj_k = eval_local_du_fluxproj(j,q_,k);

        term1 = du_flux_jp1 * Lk_p1 - du_flux_jm1 *Lk_m1;
        term2 =  e_penalty * r_lift
                * ( u_sol_jp1 * Lk_p1 - u_sol_jm1 * Lk_m1) / hjj;
        term3 = - du_proj_k ;
        term4 = - 0.5 * fact_ * ( u_sol_jp1 * dLk_p1 + u_sol_jm1 * dLk_m1 );

        viscous_resid =  simdata_->thermal_diffus
                          * ( term1 + term2 + term3 + term4 ) ;

        // Inviscid Residual Computing:
        f_proj_k = eval_localInviscidFlux_proj(q_,k);

        inviscid_resid =  - ( flux_jp1 * Lk_p1 - flux_jm1 *Lk_m1) + f_proj_k ;

        resid_[k] = fact_ * Mkk * ( viscous_resid + inviscid_resid ) ;
    }

    return;
}

double DGSolverAdvecDiffus::Compute_common_flux(const double &Ql, const double &Qr
                                      , const double &wave_speed
                                      , const double &upwind_Beta_){

    double f_upw=0.0, f_cent=0.0, f_common_=0.0;
    double aa=wave_speed;
    double BB = upwind_Beta_;

    f_upw = Rusanov(Ql,Qr);

    // f_upw = 0.5 * ( aa*(Ql+Qr) - fabs(aa) * (Qr-Ql) );

    double Fl=0.0,Fr=0.0;

    Fl = 0.5 *  pow(Ql,2);
    Fr = 0.5 *  pow(Qr,2);

    f_cent = 0.5 * (Fl+Fr) ;

    f_common_ = (BB * f_upw) + ((1.0-BB) * f_cent );

    return f_common_;
}

double DGSolverAdvecDiffus::Rusanov(const double &Ql, const double &Qr){

    int i,j;

    double Lambda_max = 0.0;
    double Fl=0.0,Fr=0.0;

    Fl = 0.5 *  pow(Ql,2);
    Fr = 0.5 *  pow(Qr,2);
    Lambda_max = 0.5 * fabs(Ql+Qr);


    return  ( 0.5 * (Fl+Fr) - 0.5 * Lambda_max * (Qr-Ql) );
}


double DGSolverAdvecDiffus::eval_burgers_inviscidFlux(const double& xi_pt_, const double *q_){

    int k;

    double xx=xi_pt_;
    double Q_;

    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_basis_poly(xx,k);

    return ( 0.5 * pow(Q_,2) );
}


double DGSolverAdvecDiffus::eval_localInviscidFlux_proj(const double *q_, const int &basis_k_){

    int k=basis_k_;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    int i;

    double II=0.0, dL_k=0.0;
    double burger_flux_=0.0;


    II=0.0;

    for (i=0; i<quad_.Nq; i++){

        dL_k = eval_basis_poly_derivative(quad_.Gaus_pts[i],k);

        burger_flux_ = eval_burgers_inviscidFlux(quad_.Gaus_pts[i],q_);

        II += quad_.Gaus_wts[i] * burger_flux_ * dL_k;
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

double DGSolverAdvecDiffus::Compute_common_sol_jump(const double &ul_
                                               , const double &ur_){
    return (ur_ - ul_);
}

double DGSolverAdvecDiffus::Compute_common_du_flux(const double& dul_
                                                  , const double& dur_){

    if(simdata_->diffus_scheme_type_=="BR2"
            || simdata_->diffus_scheme_type_=="SIP"){

        return ( 0.5 * ( dul_ + dur_) );

    }else{
        _notImplemented("diffusion Scheme type");
        return 0.0;
    }
}

void DGSolverAdvecDiffus::Compute_vertex_sol(){

    register int i;

    double Ql=0.0, Qr=0.0;

    // left boundary value:
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

void DGSolverAdvecDiffus::Compute_cont_sol(){




    return;
}

void DGSolverAdvecDiffus::Compute_exact_vertex_sol(){

    register int j;

    double xx=0.0;
    double x0,x1;

    for(j=0; j<grid_->N_exact_ppts; j++){

        xx = grid_->x_exact_ppts[j]- exact_sol_shift;

        if(simdata_->wave_form_==0){
            Q_exact[j] = eval_init_sol(xx);

        }else if(simdata_->wave_form_==1){

            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);

            if(x0==0 && x1==0)
                Q_exact[j] = 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
            else
                Q_exact[j] = (eval_init_sol(x0)+ eval_init_sol(x1));
        }
    }

    return;
}

void DGSolverAdvecDiffus::Compute_projected_exact_sol(){

    register int j; int k=0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qex_proj[j][k] = ExactSol_legendre_proj(j,k,quad_);

    return;
}

double DGSolverAdvecDiffus::eval_init_sol(const double& xx){

    double L = fabs(wave_length_);

    if(simdata_->wave_form_==0){

        return sin( simdata_->wave_freq_*PI*xx / L);

        //return ( -2. * tanh(9*xx) );

    }else if(simdata_->wave_form_==1){

        return exp(-simdata_->Gaussian_exponent_*pow(xx,2));
    }else{
        _notImplemented("Wave form is not implemented");
    }
}

double DGSolverAdvecDiffus::eval_exact_sol(double &xx){

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
    }else{
        _notImplemented("Wave form is not implemented");
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
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"polynomial order of: %d",basis_k_);
        _notImplemented(ss);
        break;
    }

//    if(basis_k_==0) {
//        return 1.0;

//    }else if(basis_k_==1) {
//        return xi_;

//    }else if(basis_k_==2) {
//        return 0.5 * (3.0 * xi_*xi_ -1.0);

//    }else if(basis_k_==3) {
//        return 0.5 * (5.0 * pow(xi_,3) - 3.0 * xi_);

//    }else if(basis_k_==4) {
//        return  (35 * pow(xi_,4) - 30 * xi_*xi_ + 3 ) /8. ;

//    }else if(basis_k_==5) {
//        return  (63.0 * pow(xi_,5) - 70.0 * pow(xi_,3) + 15.0 * xi_ ) /8. ;

//    }else {
//        char *ss=nullptr; ss= new char[100];
//        sprintf(ss,"polynomial order of: %d",basis_k_);
//        _notImplemented(ss);

//        return 0.0;
//    }

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
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"polynomial order of: %d",basis_k_);
        _notImplemented(ss);
        break;
    }

//    if(basis_k_==0) {
//        return 0.0;

//    }else if(basis_k_==1) {
//        return 1.0;

//    }else if(basis_k_==2) {
//        return (3.0 * xi_);

//    }else if(basis_k_==3) {
//        return ( 0.5 * (15.0 * pow(xi_,2) - 3.0 ) );

//    }else if(basis_k_==4) {
//        return ( 0.5 * (35.0 * pow(xi_,3) - 15.0 * xi_ ) ) ;

//    }else if(basis_k_==5) {
//        return  (315.0 * pow(xi_,4) - 210.0 * pow(xi_,2) + 15.0  ) / 8. ;

//    }else {
//        char *ss=nullptr; ss= new char[100];
//        sprintf(ss,"polynomial order of: %d",basis_k_);
//        _notImplemented(ss);

//        return 1000.0;
//    }

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

    for(k=0; k<Ndof; k++){

        Q_ += q_[k] * eval_basis_poly(xx,k);
    }

    return Q_;
}

double DGSolverAdvecDiffus::eval_local_du_fluxproj(const int eID, const double *q_
                                           , const int &basis_k_){

    int k=basis_k_;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    int i;

    double II=0.0, dL_k=0.0;
    double du_flux_=0.0;

    II=0.0;

    for (i=0; i<quad_.Nq; i++){

        dL_k = eval_basis_poly_derivative(quad_.Gaus_pts[i],k);

        du_flux_ = eval_local_du(eID, q_, quad_.Gaus_pts[i]);

        II += quad_.Gaus_wts[i] * du_flux_ * dL_k;
    }

    return II;
}

double DGSolverAdvecDiffus::ComputePolyError(){

    register int j; int i;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

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

    GaussQuad quad_;

    quad_.setup_quadrature(Ndof);

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

    GaussQuad quad_;

    quad_.setup_quadrature(Ndof);

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

    GaussQuad quad_;

    quad_.setup_quadrature(Ndof);

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

    GaussQuad quad_;

    quad_.setup_quadrature(Ndof);

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

    GaussQuad quad_;

    quad_.setup_quadrature(5);

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

    GaussQuad quad_;

    quad_.setup_quadrature(5);

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

double DGSolverAdvecDiffus::L2_error_nodal_disc_sol(){

    // We have two options: either use only
    // the interface nodes or also use the center node

    double L2_error=0.0;

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

    if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%snodal/u_cont_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,e_penalty
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
                ,e_penalty
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nfaces; j++)
            fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qv[j]);

        fclose(sol_out);
        emptyarray(fname);
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

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%saver/u_aver_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,e_penalty
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
                ,e_penalty
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
                ,e_penalty
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
                ,e_penalty
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
                ,e_penalty
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
                 ,e_penalty
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
                ,e_penalty
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
                ,e_penalty
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

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%snodal/u_disc_N%d_dt%1.3e_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,e_penalty
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
                ,e_penalty
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

    register int j; int k; int eID_=0;

    double xx=0.0,xxi_=0.0,qq=0.0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    char *fname=nullptr;
    fname = new char[200];

    // Dump time accurate continuous equally spaced solution data:
    sprintf(fname,"%stime_data/u_cont_N%d_CFL%1.3e_Beta%1.2f_Eps%1.2f_%1.3ft.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,CFL
            ,simdata_->upwind_param_
            ,e_penalty
            ,phy_time);

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->N_equal_spaced+1; j++){

        qq=0.0;

        eID_ = grid_->x_dump_to_elem[j];

        if(eID_>=0){
            xxi_ = 2.0 * (grid_->X_dump[j] - grid_->Xc[eID_])
                    / grid_->h_j[eID_];
            qq = evalSolution(&Qn[eID_][0],xxi_);
        }else{
            eID_=grid_->x_dump_to_elem[j-1];
            xxi_ = 2 * (grid_->X_dump[j] - grid_->Xc[eID_])
                    / grid_->h_j[eID_];
            qq = evalSolution(&Qn[eID_][0],xxi_);

            eID_=grid_->x_dump_to_elem[j+1];
            xxi_ = 2 * (grid_->X_dump[j] - grid_->Xc[eID_])
                    / grid_->h_j[eID_];
            qq += evalSolution(&Qn[eID_][0],xxi_);
            qq = 0.5*qq;
        }

        fprintf(sol_out,"%2.10e %2.10e\n",grid_->X_dump[j],qq);
    }

    fclose(sol_out);
    emptyarray(fname);

    // Dump time accurate Discontinuous data:
    fname = new char[200];
    sprintf(fname,"%stime_data/u_disc_N%d_CFL%1.3e_Beta%1.2f_Eps%1.3f_%1.3ft.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,CFL
            ,simdata_->upwind_param_
            ,e_penalty
            ,phy_time);

    sol_out=fopen(fname,"w");

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



    return;
}
























