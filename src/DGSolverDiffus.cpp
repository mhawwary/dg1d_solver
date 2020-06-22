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
        Nquad_viscFlux_=3;
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
    u_sol_jump = new double[grid_->Nfaces+2]; // we added 2 to account for periodic B.C. and non compactness of BR1
    Q_cont_sol = new double[grid_->N_uniform_pts];

    eta_penalty = simdata_->penalty_param_; // penalty param
    C_lift = pow(Ndof,2)/2.0; // C lifting term multiplier for LDG, BR1, BR2/SIP
    if(simdata_->diffus_scheme_type_=="SIP"
            || simdata_->diffus_scheme_type_=="BR2"){
        C_lift = eta_penalty * C_lift;
        C1_lift=0.0;
    }else if(simdata_->diffus_scheme_type_=="BR1"){
        C_lift = (1+eta_penalty)*C_lift;
        C1_lift = pow(-1,Ndof-1)*Ndof/4.0; // C1 lifting term multiplier for BR1
    }else if(simdata_->diffus_scheme_type_=="LDG"){
        C_lift = 2.0 * C_lift + eta_penalty;
        C1_lift=0.0;
    }else if(simdata_->diffus_scheme_type_=="CGR"){
        load_cgr_scheme_matrices();
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
    emptyarray(Q_cont_sol);
    emptyarray(Ndof,Lk);
    emptyarray(Ndof,dLk);
    emptyarray(Lk_norm_squar);

    quad_.Reset_quad();
    quad_viscF_.Reset_quad();

    grid_->Reset_();
    simdata_->Reset();

    return;
}

void DGSolverDiffus::load_cgr_scheme_matrices(){

    _print_log("loading CGR matrices");

    std::string scheme_name_;
    std::ostringstream oss;
    oss<<"CGR_chi"<<std::fixed<<std::setprecision(2)<<eta_penalty
          <<"_p"<<Ndof-1;
    scheme_name_ = oss.str();

    printf("cgr data fname:%s\n",scheme_name_.c_str());

    std::string fname_;
    fname_ = "./schemes_data/CGR_michigan/"+scheme_name_+".bn";
    std::ifstream sol_reader_;
    sol_reader_.open(fname_.c_str(),std::ios::in | std::ios::binary);

    // first: reading the matrix size & allocating
    int nrow_o,ncol_o;
    std::vector<std::vector<long double>> Ktemp_;
    sol_reader_.read((char*) &nrow_o,sizeof(int));
    sol_reader_.read((char*) &ncol_o,sizeof(int));
    Ktemp_.resize(nrow_o);
    Km1_s.resize(nrow_o);
    K0_s.resize(nrow_o);
    Kp1_s.resize(nrow_o);

    // second: reading the matrices by row:
    for(int i=0; i<nrow_o; i++){
        Ktemp_[i].resize(ncol_o);
        Km1_s[i].resize(ncol_o);
        sol_reader_.read((char*) &Ktemp_[i][0], ncol_o*sizeof(long double));
        for(int j=0; j<ncol_o; j++)
            Km1_s[i][j]=static_cast<double>(Ktemp_[i][j]);
    }
    for(int i=0; i<nrow_o; i++){
        Ktemp_[i].resize(ncol_o);
        K0_s[i].resize(ncol_o);
        sol_reader_.read((char*) &Ktemp_[i][0], ncol_o*sizeof(long double));
        for(int j=0; j<ncol_o; j++)
            K0_s[i][j]=static_cast<double>(Ktemp_[i][j]);
    }
    for(int i=0; i<nrow_o; i++){
        Ktemp_[i].resize(ncol_o);
        Kp1_s[i].resize(ncol_o);
        sol_reader_.read((char*) &Ktemp_[i][0], ncol_o*sizeof(long double));
        for(int j=0; j<ncol_o; j++)
            Kp1_s[i][j]=static_cast<double>(Ktemp_[i][j]);
    }
    sol_reader_.close();


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
    cout <<"===============================================\n";
    //cout << "max eigenvalue : "<<max_eigen_advec<<endl;
    //cout << "TotalVariation : "<<TV_<<endl;
    cout << "Wave length    : "<< wave_length_<<endl;
    cout << "ThermalDiffusiv: "<< simdata_->thermal_diffus << endl;
    cout << "CFL no.        : "<<CFL<<endl;
    cout << "Time step, dt  : "<<time_step<<endl;
    cout << "Last_time_step : "<<last_time_step<<endl;
    cout << "Input Nperiods : "<<simdata_->Nperiods<<endl;
    cout << "New   Nperiods : "<<simdata_->t_end_/T_period<<endl;
    cout << "Exact_sol_shift: "<<exact_sol_shift<<endl;
    cout << "T_period       : "<<T_period<<endl;
    printf( "Actual_end_time:%1.5f\n",simdata_->t_end_);
    cout << "Max_iter       : "<<simdata_->maxIter_<<endl;

    cout << "\nNumber of Elements: "<< grid_->Nelem<<"  dx:  "<<grid_->dx<<endl;
    cout << "Polynomial  order   : "<< simdata_->poly_order_  << endl;
    cout << "Runge-Kutta order   : "<< simdata_->RK_order_    << endl;
    cout << "Poly GaussQuad order: "<< Nquad_ << endl;
    cout << "Flux GaussQuad order: "<< Nquad_viscFlux_ << endl;
    cout << "Penalty parameter   : "<< eta_penalty << endl;
    cout << "Viscous Flux scheme : "<< simdata_->diffus_scheme_type_<<endl;
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
    init_wave_E_ = Compute_waveEnergy(Qn);  // initial wave energy of the projected solution

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
    double Qexact_=0.0;
    double Lk_=1.0;

    for (i=0; i<quad_.Nq; i++){
        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
        Qexact_ = eval_exact_sol(xx);
        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);
        II += quad_.Gaus_wts[i] * Qexact_ * Lk_ ;
    }
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverDiffus::UpdateResid(double **Resid_, double **Qn_){

    if(simdata_->diffus_scheme_type_=="CGR"){
        // left most cell, periodic B.C.:
        UpdateResidOneCell_cgr(0,&Qn_[grid_->Nelem-1][0]
                              ,&Qn_[0][0],&Qn_[1][0],&Resid_[0][0]);
        // right most cell, periodic B.C.:
        UpdateResidOneCell_cgr(grid_->Nelem-1,&Qn_[grid_->Nelem-2][0]
                              ,&Qn_[grid_->Nelem-1][0],&Qn_[0][0],&Resid_[grid_->Nelem-1][0]);
        // interior cells:
        for(register int j=1; j<grid_->Nelem-1; j++)
            UpdateResidOneCell_cgr(j,&Qn_[j-1][0],&Qn_[j][0],&Qn_[j+1][0],&Resid_[j][0]);
        return;
    }else{

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
        flux_com[j] = Compute_common_du_flux(dQl,dQr);
        u_sol_jump[j+Nghost_l] = Compute_common_sol_jump(Ql,Qr);
        flux_com[grid_->Nfaces-1] = flux_com[j];
        u_sol_jump[grid_->Nfaces-1+Nghost_l] = u_sol_jump[j+Nghost_l];

        // Interior Faces
        for(j=1; j<grid_->Nfaces-1; j++){
            Ql = evalSolution(&Qn_[j-1][0], 1);
            Qr = evalSolution(&Qn_[j][0], 0);
            dQl = eval_local_du_fast(j-1, &Qn_[j-1][0], 1);
            dQr = eval_local_du_fast(j, &Qn_[j][0], 0);
            // Viscous Common Flux:
            flux_com[j] = Compute_common_du_flux(dQl,dQr);
            u_sol_jump[j+Nghost_l] = Compute_common_sol_jump(Ql,Qr);
        }
        // Updating the additional jump values for BR1:
        u_sol_jump[0] = u_sol_jump[grid_->Nfaces-2+Nghost_l]; // for left Boundary
        u_sol_jump[grid_->Nfaces+Nghost_l] = u_sol_jump[1+Nghost_l]; // for right Boundary
        // Element loop to calculate and update the residual:
        //----------------------------------------------------
        for(j=0; j<grid_->Nelem; j++){
            UpdateResidOneCell(j, &Qn_[j][0], &Resid_[j][0]);
        }

        return;
    }
}

void DGSolverDiffus::UpdateResidOneCell(const int &cellid,double *q_, double *resid_){

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
        term2 = tempC_lift*( (u_jump_jp12 * Lk_p1) - (u_jump_jm12 * Lk_m1)); // local lift term
        term3 = tempC1_lift*( (u_jump_jp32 + u_jump_jm12 )* Lk_p1
                              -(u_jump_jp12 + u_jump_jm32 )* Lk_m1); // for BR1, 0.0 for others
        term4 = - du_proj_k ;
        if(simdata_->diffus_scheme_type_=="LDG")
            term5 = - fact_ * u_jump_jp12 * dLk_p1 ;
        else
            term5 = - 0.5 * fact_ * ( (u_jump_jp12 * dLk_p1) + (u_jump_jm12 * dLk_m1) );

        resid_[k] = fact_ * Mkk * simdata_->thermal_diffus * ( term1 + term2 + term3 + term4 + term5) ;
    }

    return;
}


void DGSolverDiffus::UpdateResidOneCell_cgr(const int &eID_,const double *qm1_,const double *q0_
                                            ,const double *qp1_,double *resid_){

    double h0_inv,hm1_inv,hp1_inv;
    if(eID_==0) //periodic B.C.
        hm1_inv=1./(grid_->h_j[grid_->Nelem-1]*grid_->h_j[grid_->Nelem-1]);
    else
        hm1_inv = 1./(grid_->h_j[eID_-1]*grid_->h_j[eID_-1]);

    h0_inv  = 1./(grid_->h_j[eID_]*grid_->h_j[eID_]);

    if(eID_==grid_->Nelem-1) //periodic B.C.
        hp1_inv=1./(grid_->h_j[0]*grid_->h_j[0]);
    else
        hp1_inv = 1./(grid_->h_j[eID_+1]*grid_->h_j[eID_+1]);


    double temp_m1,temp_0,temp_p1;
    for(int k=0; k<Ndof; k++){
        temp_m1=0.; temp_0=0.; temp_p1=0.;
        for(int j=0; j<Ndof; j++){
            temp_m1+= Km1_s[k][j]*qm1_[j];
            temp_0 += K0_s [k][j]*q0_[j];
            temp_p1+= Kp1_s[k][j]*qp1_[j];
        }
        resid_[k]=simdata_->thermal_diffus*(temp_m1*hm1_inv,temp_0*h0_inv+temp_p1*hp1_inv);
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

double DGSolverDiffus::Compute_common_sol_jump(const double ul_
                                               , const double ur_){
    double jump_ = ur_ - ul_;
    return jump_;
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

    for(j=0; j<grid_->N_exact_ppts; j++){
        xx = grid_->x_exact_ppts[j];
        Q_exact[j] = eval_exact_sol(xx);
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
    double wave_value_=0.0;

    if(simdata_->wave_form_==0){  // u(x,0) = A * sin/cos ( (f * PI * x) / L + phy ) + C
        double argum_ = (simdata_->wave_freq_*PI*xx)
                / fabs(wave_length_) + simdata_->wave_shift ;
        if(simdata_->trig_wave_type_=="sin")
            wave_value_ = simdata_->wave_amp_ * sin(argum_)
                    + simdata_->wave_const;
        else if(simdata_->trig_wave_type_=="cos")
            wave_value_ = simdata_->wave_amp_ * cos(argum_)
                    + simdata_->wave_const;

    }else if(simdata_->wave_form_==1){   // Gaussian wave, f(x) = Amp * exp(-b*x^2)
        wave_value_ = simdata_->Gaussian_amp_
                * exp(-simdata_->Gaussian_exponent_*pow(xx,2)) ;
    }else{
        _notImplemented("Wave form is not implemented");
    }
    return wave_value_;
}

double DGSolverDiffus::eval_exact_sol(double &xx_){

    double Qexx_=0.0;

    if(simdata_->wave_form_==0){  // single wave mode either a sine or a cosine
        double wave_initial_value = eval_init_sol(xx_);
        double Kfreq_ = pow(simdata_->wave_freq_*PI/fabs(wave_length_), 2);
        Qexx_ = wave_initial_value
                * exp(-Kfreq_ * simdata_->thermal_diffus *phy_time);
    }else if(simdata_->wave_form_==1){ // Gaussian wave, f(x) = Amp * exp(-b*x^2)
        register int j;
        double b_ = simdata_->Gaussian_exponent_;
        double gamma_=simdata_->thermal_diffus;
        double L_ = 0.5*wave_length_;  // half of periodic domain width
        double sqrtb_L_ = sqrt(b_) * L_;
        double sqrtPib_ = sqrt(PI/b_);
        static int n_eigenfunc_gaussian_ = 500; // initial number of eigenvalues to use for the exact solution expansion
        int n_eignfunc_temp=n_eigenfunc_gaussian_;

        double A0_,A_,m;
        std::complex<double> z_;

        A0_ =(sqrtPib_/(2.0*L_)) * Faddeeva::erf(sqrtb_L_); // Faddeeva is used for erf(complex numbers)
        Qexx_ = A0_;
        for(j=1; j<n_eigenfunc_gaussian_+1; j++){
            m = j*PI/L_;
            z_ = sqrtb_L_ + 1i*j*PI/(2.0*sqrtb_L_);
            A_ = (sqrtPib_/L_) * exp(-m*m/(4.0*b_)) * real(Faddeeva::erf(z_));
            Qexx_ += A_ * cos(m*xx_) * exp(-gamma_* m*m *phy_time);
            if(A_ <= 1e-10){ // just some tolerance
                n_eigenfunc_gaussian_=j;
                j=1e6;
                break;
            }
        }
        Qexx_ = Qexx_ * simdata_->Gaussian_amp_;

        if(n_eigenfunc_gaussian_ != n_eignfunc_temp)
            cout << "\nNumber of eigen functions for the Gaussian exact solution: "
                 << n_eigenfunc_gaussian_ <<endl;
    }else{
        _notImplemented("Wave form is not implemented");
    }

    return Qexx_;
}

double DGSolverDiffus::eval_basis_poly(const double& xi_, const int& basis_k_){

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
    double L1_error=0.0,elem_error=0.0,II=0.0,q_ex,q_n;

    for(j=0; j<grid_->Nelem; j++){
        elem_error=0.0;
        for(i=0; i<quad_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j][0], quad_.Gaus_pts[i]);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            elem_error += quad_.Gaus_wts[i] * fabs(q_ex - q_n);
        }
        II += ( grid_->h_j[j] * elem_error) ;
    }
    L1_error = 0.5 *II/(grid_->xf-grid_->x0);

    return L1_error;
}

double DGSolverDiffus::L2_error_projected_sol(){

    register int j; int i;
    double L2_error=0.0,elem_error=0.0,II=0.0,q_ex,q_n;
    double xx=0.0;

    for(j=0; j<grid_->Nelem; j++){
        elem_error=0.0;
        for(i=0; i<quad_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j][0], quad_.Gaus_pts[i]);
            //xx =  0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
            //q_ex = eval_exact_sol(xx);
            q_n = evalSolution(&Qn[j][0],quad_.Gaus_pts[i]);
            elem_error += quad_.Gaus_wts[i] * pow((q_ex - q_n),2);
        }
        II += (  grid_->h_j[j] * elem_error) ;
        //printf("eID: %d, \t err: %1.5f",j,elem_error);
        //std::cin.get();
    }
    L2_error = sqrt(0.5 *II/(grid_->xf-grid_->x0));

    return L2_error;
}

double DGSolverDiffus::L1_error_average_sol(){

    register int j;
    double L1_error=0.0,error=0.0;
    for(j=0; j<grid_->Nelem; j++)
        error += fabs(Qex_proj[j][0] - Qn[j][0]);
    L1_error = error/grid_->Nelem;

    return L1_error;
}

double DGSolverDiffus::L2_error_average_sol(){

    register int j;
    double L2_error=0.0,error=0.0;
    for(j=0; j<grid_->Nelem; j++)
        error += pow((Qex_proj[j][0] - Qn[j][0]),2);
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
    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"%serrors/errors_N%d_CFL%1.4f_Eps%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,CFL
            ,eta_penalty
            ,simdata_->Nperiods);

    FILE* solerror_out=nullptr;

    if(phy_time==0 && simdata_->restart_flag==1)
        solerror_out=fopen(fname,"at+");
    else if(phy_time==0 && simdata_->restart_flag==0)
        solerror_out=fopen(fname,"w");
    else
        solerror_out=fopen(fname,"at+");

    double wave_energy_ = Compute_waveEnergy(Qn);
    double wave_exact_energy_ = Compute_waveEnergy(Qex_proj);
    double GG_ = wave_energy_/init_wave_E_;
    double GG_ex_ = wave_exact_energy_/init_wave_E_;
    double L1_proj_sol_ = L1_error_projected_sol();
    double L2_proj_sol_ = L2_error_projected_sol();
    //double L2_proj_sol_ = L2_error_modal();
    double L1_nodal_sol_ = L1_error_nodal_gausspts_proj();
    double L2_nodal_sol_ = L2_error_nodal_gausspts_proj();
    //double L1_nodal_sol_ = L1_error_nodal_cont_sol(); //for testing computing at the same nodes as FD
    //double L2_nodal_sol_ = L2_error_nodal_cont_sol(); // but need to make N_exact = N_uniform in griddata
    double L1_aver_sol_ = L1_error_average_sol();
    double L2_aver_sol_ = L2_error_average_sol();

    fprintf(solerror_out, "%1.10f %2.10e %2.10e %2.10e %2.10e\n"
            ,phy_time,L1_proj_sol_, L2_proj_sol_
            , L1_aver_sol_, L2_aver_sol_);

    fclose(solerror_out);
    emptyarray(fname);

    printf("  L1_proj_sol:%2.5e  L2_proj_sol:%2.5e  L1_aver:%2.5e L2_aver:%2.5e E_ex:%1.5f E_num:%1.5f G_ex:%1.5f G_num:%1.5f"
           ,L1_proj_sol_, L2_proj_sol_, L1_aver_sol_, L2_aver_sol_
           ,wave_exact_energy_,wave_energy_,GG_ex_,GG_);
    dump_timeaccurate_waveenergy(wave_exact_energy_,wave_energy_,GG_ex_,GG_);

    return;
}

void DGSolverDiffus::dump_timeaccurate_waveenergy(const double& in_E_ex_,
                         const double& in_E_,
                         const double& in_GG_ex_,
                         const double& in_GG_){
    char *fname=nullptr;
    fname = new char[100];
    if(simdata_->Sim_mode=="CFL_const"
            || simdata_->Sim_mode =="error_analysis_CFL"
            || simdata_->Sim_mode =="test"
            || simdata_->Sim_mode =="normal")
        sprintf(fname,"%stime_data/wave_energy_N%d_CFL%1.4f_Eps%1.2f.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty);
    else if(simdata_->Sim_mode=="dt_const"
            || simdata_->Sim_mode=="error_analysis_dt" )
        sprintf(fname,"%stime_data/wave_energy_N%d_dt%1.3e_Eps%1.2f.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,eta_penalty);

    FILE* sol_out=nullptr;

    if(phy_time==0 && simdata_->restart_flag==1)
        sol_out=fopen(fname,"at+");
    else if(phy_time==0 && simdata_->restart_flag==0)
        sol_out=fopen(fname,"w");
    else
        sol_out=fopen(fname,"at+");

    fprintf(sol_out, "%1.10f %1.10e %1.10e %1.5f %1.5f\n"
            ,phy_time, in_E_ex_, in_E_, in_GG_ex_, in_GG_);

    fclose(sol_out);
    emptyarray(fname);
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
        sprintf(fname,"%stime_data/u_cont_N%d_CFL%1.4f_Eps%1.2f_%1.4ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        // Dump time accurate continuous equally spaced solution data:
        sprintf(fname,"%stime_data/u_cont_N%d_dt%1.3e_Eps%1.2f_%1.4ft.dat"
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
        sprintf(fname,"%stime_data/u_disc_N%d_CFL%1.4f_Eps%1.2f_%1.4ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        sprintf(fname,"%stime_data/u_disc_N%d_dt%1.3e_Eps%1.2f_%1.4ft.dat"
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
    //===========================================
    Compute_projected_exact_sol();
    Compute_exact_vertex_sol();

    //Discontinuous exact solution:
    fname = new char[100];
    sprintf(fname,"%stime_data/u_disc_exact_N%d_%1.4ft.dat"
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

    ////Continuous Exact solution:
    //Compute_TimeAccurate_exact_sol();
    fname = new char[100];
    sprintf(fname,"%stime_data/u_cont_exact_%1.4ft.dat"
            ,simdata_->case_postproc_dir
            ,phy_time);

    sol_out=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact[j]);

    fclose(sol_out);
    emptyarray(fname);

    //dump_gaussian_coeffs(); // was a trial to dump exact wavenumbers amplitudes

    return;
}

void DGSolverDiffus::compute_uniform_cont_sol(){

    //_print("-------Inside compute unifrom ");

    // Dump continuous data on uniform points:
    //-------------------------------------------
    register int j; int k;

    int count_=0; //continuous points counter

    // For first grid point, the first face:
    //----------------------------------------
    k = simdata_->N_uniform_pts_per_elem_-1;
    j= grid_->Nelem-1;   // solution from left cell
    Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]); // a potential bug, fixed by changing j-1 to j

    k=0; j=0; //solution from right cell
    Q_cont_sol[count_] += evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
    Q_cont_sol[count_] = 0.5*Q_cont_sol[count_]; //averaging
    count_++;

    // Element zero interior points:
    //-------------------------------
    for(k=1; k<simdata_->N_uniform_pts_per_elem_-1; k++) {
        Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
        count_++;
    }

    // For the rest of the faces/points:
    //-----------------------------------
    for(j=1; j<grid_->Nelem; j++){

        //left interface point for elment j:
        k=0;
        Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
        k = simdata_->N_uniform_pts_per_elem_-1;
        Q_cont_sol[count_] += evalSolution(&Qn[j-1][0],grid_->xi_uniform[k]);
        Q_cont_sol[count_] = 0.5*Q_cont_sol[count_]; // averaging
        count_++;

        //interior points for elment j:
        for(k=1; k<simdata_->N_uniform_pts_per_elem_-1; k++) {
            Q_cont_sol[count_] = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
            count_++;
        }
    }
    Q_cont_sol[count_] = Q_cont_sol[0]; //final right interface/boundary, periodic B.C.

    //printf("\nj:%d \t Count: %d \t N_uniform: %d\n",j, count_, grid_->N_uniform_pts);
    //std::cin.get();
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

double DGSolverDiffus::Compute_waveEnergy(double **in_Qn_){
    // This is not exactly an energy quantity, it is precisely:
    // an averaged amplitude.
    // So, E = (wave_energy)^2 and KE = 0.5*E
    register int j; int i;
    double wave_energy=0.0,elem_energy=0.0,II=0.0,q_n;

    for(j=0; j<grid_->Nelem; j++){
        elem_energy=0.0;
        for(i=0; i<quad_.Nq; i++) {
            q_n = evalSolution(&in_Qn_[j][0],quad_.Gaus_pts[i]);
            elem_energy += quad_.Gaus_wts[i] * q_n*q_n;
        }
        II +=   grid_->h_j[j] * elem_energy ;
    }
    wave_energy = sqrt(0.5*II/(grid_->xf-grid_->x0));

    return wave_energy;
}

double DGSolverDiffus::L1_norm_modal(double **Qn_, const int in_Nelem_
                                     , const int in_ndof_){

    register int j; int k;
    double L1_norm=0.0,II=0.0;

    for(j=0; j<in_Nelem_; j++){
        for(k=0; k<in_ndof_; k++) {
            II += fabs(Qn_[j][k]);
        }
    }
    L1_norm = II/(in_ndof_*in_Nelem_);

    return L1_norm;
}

double DGSolverDiffus::L2_norm_modal(double **Qn_, const int in_Nelem_
                                     , const int in_ndof_){
    register int j; int k;
    double L2_norm=0.0,II=0.0;

    for(j=0; j<in_Nelem_; j++){
        for(k=0; k<in_ndof_; k++) {
            II += pow(Qn_[j][k],2);
        }
    }
    L2_norm = sqrt(II/(in_ndof_*in_Nelem_));

    return L2_norm;
}

double DGSolverDiffus::L2_error_modal(){
    register int j; int k;
    double L2_error=0.0,II=0.0;

    for(j=0; j<grid_->Nelem; j++){
        for(k=0; k<Ndof; k++) {
            II += pow(Qn[j][k]-Qex_proj[j][k],2);
        }
    }
    L2_error = sqrt(II/(Ndof*grid_->Nelem));

    return L2_error;
}

void DGSolverDiffus::dump_gaussian_coeffs(){

    char *fname=nullptr;
    fname = new char[250];

    if(simdata_->Sim_mode=="CFL_const"
            || simdata_->Sim_mode =="error_analysis_CFL"
            || simdata_->Sim_mode =="test"
            || simdata_->Sim_mode =="normal"){
        sprintf(fname,"%stime_data/gaussian_amp_N%d_CFL%1.4f_Eps%1.2f_%1.4ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,eta_penalty
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        sprintf(fname,"%stime_data/gaussian_amp_N%d_dt%1.3e_Eps%1.2f_%1.4ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,eta_penalty
                ,phy_time);
    }

    FILE* sol_out=fopen(fname,"w");

    register int j;
    double b_ = simdata_->Gaussian_exponent_;
    double gamma_=simdata_->thermal_diffus;
    double L_ = 0.5*wave_length_;  // half of periodic domain width
    double sqrtb_L_ = sqrt(b_) * L_;
    double sqrtPib_ = sqrt(PI/b_);
    static int n_eigenfunc_gaussian_ = 500; // initial number of eigenvalues to use for the exact solution expansion
    int n_eignfunc_temp=n_eigenfunc_gaussian_;

    double A0_,A_,m,Amp_;
    std::complex<double> z_;

    A0_ =(sqrtPib_/(2.0*L_)) * Faddeeva::erf(sqrtb_L_); // Faddeeva is used for erf(complex numbers)
    fprintf(sol_out,"%d %2.10e\n",0,A0_);

    for(j=1; j<n_eigenfunc_gaussian_+1; j++){
        m = j*PI/L_;
        z_ = sqrtb_L_ + 1i*j*PI/(2.0*sqrtb_L_);
        A_ = (sqrtPib_/L_) * exp(-m*m/(4.0*b_)) * real(Faddeeva::erf(z_));
        Amp_= A_ * exp(-gamma_* m*m *phy_time);

        fprintf(sol_out,"%d %2.10e\n",j,Amp_);

        if(A_ <= 1e-10){ // just some tolerance
            n_eigenfunc_gaussian_=j;
            j=1e6;
            break;
        }
    }

    // error checking message
    if(n_eigenfunc_gaussian_ != n_eignfunc_temp)
        cout << "\nNumber of eigen functions for the Gaussian exact solution: "
             << n_eigenfunc_gaussian_ <<endl;

    fclose(sol_out);
    emptyarray(fname);

    return;
}
