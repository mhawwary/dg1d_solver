#include "DGSolverDiffus.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

DGSolverDiffus::~DGSolverDiffus(void){

    Reset_solver();
}

void DGSolverDiffus::setup_solver(GridData& meshdata_, SimData& osimdata_){

    register int i;

    grid_ = &meshdata_;
    simdata_ = &osimdata_;

    Ndof = simdata_->poly_order_+1;

    // Nquad  is for error integrations and initialization
    // Nquad_invFlux & Nquad_viscFlux_ are for flux projection
    switch (Ndof) {
    case 1:  // p0
        Nquad_ = 8;
        Nquad_viscFlux_=1; // it is not applicable
        break;
    case 2: // p1
        Nquad_ = 8;
        Nquad_viscFlux_=1;
        break;
    case 3:  // p2
        Nquad_ = 8;
        Nquad_viscFlux_=2;
        break;
    case 4:  // p3
        Nquad_ = 8;
        Nquad_viscFlux_=3;
        break;
    case 5:  // p4
        Nquad_ = 8;
        Nquad_viscFlux_=4;
        break;
    case 6:  // p5
        Nquad_ = 8;
        Nquad_viscFlux_=5;
        break;
    default:
        break;
    }

    quad_.setup_quadrature(Nquad_);
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
    u_sol_jump = new double[grid_->Nfaces];
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

void DGSolverDiffus::setup_basis_interpolation_matrices(){

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

void DGSolverDiffus::Reset_solver(){

    emptyarray(grid_->Nelem,Qn);
    emptyarray(Q_exact);
    emptyarray(flux_com);
    emptyarray(u_sol_jump);
    emptyarray(Qv);
    emptyarray(grid_->Nelem,Qex_proj);

    grid_->Reset_();


    return;
}


// Solver functions
//-------------------------------------------

void DGSolverDiffus::CalcTimeStep(){

    compute_uniform_cont_sol();
    //double TV_ = compute_totalVariation();

    double dx2 = pow(grid_->dx,2);
    double radius_diffus_ =0.0;

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;
    radius_diffus_ = simdata_->thermal_diffus / dx2 ;

    if(simdata_->calc_dt_flag==1){ // use CFL as input
        CFL = simdata_->CFL_;
        time_step = CFL / radius_diffus_;
        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){  // use dt as input
        time_step = simdata_->dt_;
        last_time_step = time_step;
        CFL = time_step * radius_diffus_ ;
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
    //cout << "max eigenvalue : "<<max_eigen_advec<<endl;
    //cout << "TotalVariation : "<<TV_<<endl;
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
    cout << "Penalty parameter : "<< eta_penalty << endl;
    cout << "Poly GaussQuad order  : "<< Nquad_ << endl;
    cout << "Flux GaussQuad order  : "<< Nquad_viscFlux_ << endl;
    cout <<"===============================================\n";

    return;
}

void DGSolverDiffus::InitSol(){

    register int j;
    int k=0;
    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qn[j][k] = initSol_legendre_proj(j,k,quad_);

    // max_eigen_diffus =  // if needed

    CalcTimeStep(); // based on maximum eigenvalues if any

    return;
}

double DGSolverDiffus::initSol_legendre_proj(const int &eID,
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
        Qinit_= eval_init_sol(xx);
        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);
        II += quad_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverDiffus::ComputeExactSolShift(){
    exact_sol_shift = (wave_speed_ * phy_time );
    return;
}

double DGSolverDiffus::ExactSol_legendre_proj(const int &eID,
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

        if(simdata_->wave_form_==0){       // single mode wave
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
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverDiffus::UpdateResid(double **Resid_, double **Qn_){

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
    flux_com[j] = Compute_common_du_flux(dQl,dQr);
    u_sol_jump[j] = Compute_common_sol_jump(Ql,Qr);
    flux_com[grid_->Nfaces-1] = flux_com[j];
    u_sol_jump[grid_->Nfaces-1] = u_sol_jump[j];

    // Interior Faces
    for(j=1; j<grid_->Nfaces-1; j++){
        Ql = evalSolution(&Qn_[j-1][0], 1);
        Qr = evalSolution(&Qn_[j][0], 0);
        dQl = eval_local_du_fast(j-1, &Qn_[j-1][0], 1);
        dQr = eval_local_du_fast(j, &Qn_[j][0], 0);
        // Viscous Common Flux:
        flux_com[j] = Compute_common_du_flux(dQl,dQr);
        u_sol_jump[j] = Compute_common_sol_jump(Ql,Qr);
    }
    // Element loop to calculate and update the residual:
    //----------------------------------------------------
    for(j=0; j<grid_->Nelem; j++){
        UpdateResidOneCell(j, &Qn_[j][0], &Resid_[j][0]);
    }

    return;
}

void DGSolverDiffus::UpdateResidOneCell(const int &cellid, double *q_, double *resid_){

    // General parameters:
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

    if(simdata_->diffus_scheme_type_=="LDG"){
        // giving that beta_e+1/2 = 1, beta_e-1/2=0
        u_jump_jm12 = 0.0;
        u_jump_jp12 = 2.0*u_sol_jump[j+1];
    }else if(simdata_->diffus_scheme_type_=="SIP"
             ||simdata_->diffus_scheme_type_=="BR2"){
        u_jump_jm12 = u_sol_jump[j];
        u_jump_jp12 = u_sol_jump[j+1];
    }else if(simdata_->diffus_scheme_type_=="BR1"){
        u_jump_jm12 = u_sol_jump[j];
        u_jump_jp12 = u_sol_jump[j+1];
        u_jump_jp32 = u_sol_jump[j+2];
        u_jump_jm32 = u_sol_jump[j-1];
    }

    du_flux_jm1 = flux_com[j];
    du_flux_jp1 = flux_com[j+1];

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
        term5 = - 0.5 * fact_ * ( u_jump_jp12 * dLk_p1 + u_jump_jm12 * dLk_m1 );

        resid_[k] = fact_ * Mkk
                * simdata_->thermal_diffus
                * ( term1 + term2 + term3 + term4 + term5) ;
    }

    return;
}

double DGSolverDiffus::eval_local_du(const int eID,
                                     const double *q_
                                     , const double& xi_pt){

    double fact_ = 2.0/grid_->h_j[eID];
    double II=0.0;
    int k;

    for(k=0; k<Ndof; k++)
        II += q_[k] * eval_basis_poly_derivative(xi_pt, k);

    return (fact_*II);
}

double DGSolverDiffus::eval_local_du_fast(const int eID,
                                          const double *q_
                                          , const int& i_pos_){
    double fact_ = 2.0/grid_->h_j[eID];
    double II=0.0;
    int k;

    for(k=0; k<Ndof; k++)
        II += q_[k] * dLk[k][i_pos_];

    return (fact_*II);
}

double DGSolverDiffus::Compute_common_sol_jump(const double &ul_
                                               , const double &ur_){
    return (ur_ - ul_);
}

double DGSolverDiffus::Compute_common_du_flux(const double& dul_
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

void DGSolverDiffus::Compute_vertex_sol(){

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

void DGSolverDiffus::Compute_exact_vertex_sol(){
    register int j;
    double xx=0.0;
    double x0,x1;

    ComputeExactSolShift(); // updating the shift for the currnet time

    for(j=0; j<grid_->N_exact_ppts; j++){
        xx = grid_->x_exact_ppts[j]- exact_sol_shift;

        if(simdata_->wave_form_==0){  // single wave mode
            Q_exact[j] = eval_init_sol(xx);
        }else if(simdata_->wave_form_==1){ // Gaussian wave
            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);
            if(x0==0 && x1==0)
                Q_exact[j] = 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
            else
                Q_exact[j] = (eval_init_sol(x0)+ eval_init_sol(x1));
        }else{
            FatalError_exit("Wave form is not implemented");
        }
    }

    return;
}

void DGSolverDiffus::Compute_projected_exact_sol(){
    register int j; int k=0;
    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qex_proj[j][k] = ExactSol_legendre_proj(j,k,quad_);
    return;
}

double DGSolverDiffus::eval_init_sol(const double& xx){

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

double DGSolverDiffus::eval_exact_sol(double &xx){

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

double DGSolverDiffus::eval_basis_poly(const double& xi_, const int& basis_k_){

    if(basis_k_==0) {
        return 1.0;

    }else if(basis_k_==1) {
        return xi_;

    }else if(basis_k_==2) {
        return 0.5 * (3.0 * xi_*xi_ -1.0);

    }else if(basis_k_==3) {
        return 0.5 * (5.0 * pow(xi_,3) - 3.0 * xi_);

    }else if(basis_k_==4) {
        return  (35 * pow(xi_,4) - 30 * xi_*xi_ + 3 ) /8. ;

    }else {
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"polynomial order of: %d",basis_k_);
        _notImplemented(ss);

        return 0.0;
    }
}

double DGSolverDiffus::eval_basis_poly_derivative(const double& xi_
                                                  , const int& basis_k_){
    /*if(basis_k_==0) {
        return 0.0;

    }else if(basis_k_==1) {
        return 1.0;

    }else if(basis_k_==2) {
        return (3.0 * xi_);

    }else if(basis_k_==3) {
        return ( 0.5 * (15.0 * pow(xi_,2) - 3.0 ) );

    }else if(basis_k_==4) {
        return ( 0.5 * (35.0 * pow(xi_,3) - 15.0 * xi_ ) ) ;

    }else {
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"polynomial order of: %d",basis_k_);
        _notImplemented(ss);

        return 0.0;
    }*/

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

double DGSolverDiffus::eval_basis_norm_squared(const int &basis_k_){

    // this is norm^2

    return 2./(2*basis_k_+1);
}

double DGSolverDiffus::evalSolution(const double* q_, const double& xi_pt_){
    double xx=xi_pt_;
    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_basis_poly(xx,k);

    return Q_;
}

double DGSolverDiffus::evalSolution(const double* q_, const int& i_pos_){

    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * Lk[k][i_pos_];   // i_pos_ = 0 if xi =-1, i_pos_ =1 if xi =1

    return Q_;
}

double DGSolverDiffus::eval_local_du_fluxproj(const int eID, const double *q_
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

double DGSolverDiffus::eval_local_du_fluxproj_exact(const int eID, const double *q_
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

double DGSolverDiffus::ComputePolyError(){

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

double DGSolverDiffus::L1_error_nodal_gausspts_proj(){

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

double DGSolverDiffus::L2_error_nodal_gausspts_proj(){

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

double DGSolverDiffus::L1_error_nodal_gausspts(){

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

double DGSolverDiffus::L2_error_nodal_gausspts(){

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

double DGSolverDiffus::L1_error_projected_sol(){

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

double DGSolverDiffus::L2_error_projected_sol(){

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

double DGSolverDiffus::L1_error_average_sol(){

    register int j;

    double L1_error=0.0,error=0.0;

    for(j=0; j<grid_->Nelem; j++){

        error += fabs(Qex_proj[j][0] - Qn[j][0]);
    }

    L1_error = error/grid_->Nelem;

    return L1_error;
}

double DGSolverDiffus::L2_error_average_sol(){

    register int j;

    double L2_error=0.0,error=0.0;

    for(j=0; j<grid_->Nelem; j++){

        error += pow((Qex_proj[j][0] - Qn[j][0]),2);
    }

    L2_error = sqrt(error/grid_->Nelem);

    return L2_error;
}

void DGSolverDiffus::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[100];

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

        sprintf(fname,"%snodal/u_cont_N%d_CFL%1.3f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
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

void DGSolverDiffus::print_average_sol(){

    register int j;

    char *fname=nullptr;
    fname = new char[100];

    GaussQuad quad_;

    quad_.setup_quadrature(5);

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

        sprintf(fname,"%saver/u_aver_N%d_CFL%1.3f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
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

void DGSolverDiffus::dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                                 ,double& L1_aver_sol_,double& L2_aver_sol_
                                 ,double& L1_nodal_gausspts, double& L2_nodal_gausspts){

    char *fname=nullptr;
    fname = new char[100];

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

void DGSolverDiffus::dump_timeaccurate_errors(){

    //Fix me! You always need to make sure that the exact solution has been updated, i.e,
    // Compute_projected_exact_sol(); Compute_exact_vertex_sol(); has been called
    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"%serrors/errors_N%d_CFL%1.4f_Eps%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,CFL
            ,eta_penalty
            ,simdata_->Nperiods);

    FILE* solerror_out=fopen(fname,"at+");

    double L1_proj_sol_ = L1_error_projected_sol();
    double L2_proj_sol_ = L2_error_projected_sol();
    double L1_nodal_sol_ = L1_error_nodal_gausspts();
    double L2_nodal_sol_ = L2_error_nodal_gausspts();
    //double L1_nodal_sol_ = L1_error_nodal_cont_sol(); //for testing computing at the same nodes as FD
    //double L2_nodal_sol_ = L2_error_nodal_cont_sol(); // but need to make N_exact = N_uniform in griddata

    fprintf(solerror_out, "%1.10f %2.10e %2.10e %2.10e %2.10e\n"
            ,phy_time,L1_proj_sol_, L2_proj_sol_
            , L1_nodal_sol_, L2_nodal_sol_);

    fclose(solerror_out);
    emptyarray(fname);

    printf("  L1_proj_sol:%2.5e  L2_proj_sol:%2.5e  L1_nodal:%2.5e  L2_nodal:%2.5e"
           ,L1_proj_sol_, L2_proj_sol_, L1_nodal_sol_, L2_nodal_sol_);

    return;
}

void DGSolverDiffus::dump_discont_sol(){

    register int j; int k;

    double xx=0.0,qq=0.0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    char *fname=nullptr;
    fname = new char[100];

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

        sprintf(fname,"%snodal/u_disc_N%d_CFL%1.3f_Eps%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
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

void DGSolverDiffus::dump_timeaccurate_sol(){

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
        sprintf(fname,"%stime_data/u_cont_N%d_CFL%1.4f_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        // Dump time accurate continuous equally spaced solution data:
        sprintf(fname,"%stime_data/u_cont_N%d_dt%1.3e_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
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
        sprintf(fname,"%stime_data/u_disc_N%d_CFL%1.4f_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        sprintf(fname,"%stime_data/u_disc_N%d_dt%1.3e_Eps%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
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

    // Dumping Exact solution data:
    Compute_projected_exact_sol();
    Compute_exact_vertex_sol();

    fname = new char[100];
    sprintf(fname,"%stime_data/u_disc_exact_N%d_%1.3ft.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,phy_time);

    sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<grid_->N_xi_disc_ppts; k++) {
            qq = evalSolution(&Qex_proj[j][0],grid_->xi_disc[k]);
            xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                    + grid_->Xc[j];

            fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
        }

    fclose(sol_out);
    emptyarray(fname);

    //==============================================
    //Compute_TimeAccurate_exact_sol();
    fname = new char[100];
    sprintf(fname,"%stime_data/u_cont_exact_%1.3ft.dat"
            ,simdata_->case_postproc_dir
            ,phy_time);

    sol_out=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact[j]);

    fclose(sol_out);
    emptyarray(fname);

    return;
}

void DGSolverDiffus::compute_uniform_cont_sol(){

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

double DGSolverDiffus::compute_totalVariation(){

    register int i;

    double TV_=0.0;
    //max_eigen_advec = 0.0;

    for(i=1; i<grid_->N_uniform_pts; i++){
        TV_ += fabs(Q_cont_sol[i]-Q_cont_sol[i-1]);
        ///if(fabs(Q_cont_sol[i-1])>max_eigen_advec)
        //  max_eigen_advec=fabs(Q_cont_sol[i-1]);
    }

    return TV_;
}
