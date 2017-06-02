#include "DGSolver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
DGSolver::DGSolver(void){}

DGSolver::~DGSolver(void){

    Reset_solver();
}

void DGSolver::setup_solver(GridData& meshdata_, SimData& osimdata_){

    grid_ = &meshdata_;
    simdata_ = &osimdata_;

    Ndof= simdata_->poly_order_+1;

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
    cout <<"===============================================\n";

    return;
}

void DGSolver::Reset_solver(){

    //emptyarray(X);
    //emptyarray(Xc);
    //emptyarray(h_j);
    emptyarray(xi);

    emptyarray(grid_->Nelem,Qn);
    emptyarray(Q_exact);
    emptyarray(flux_com);
    emptyarray(Qv);
    emptyarray(grid_->Nelem,Qex_proj);

    grid_->Reset_();


    return;
}


// Solver functions
//-------------------------------------------

void DGSolver::CalcTimeStep(){

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;

    if(simdata_->calc_dt_flag==1){

        CFL = simdata_->CFL_;
        time_step = (grid_->dx * CFL )/ simdata_->a_wave_;
        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){

        time_step = simdata_->dt_;
        last_time_step = time_step;
        CFL = simdata_->a_wave_ * time_step / grid_->dx ;
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

void DGSolver::InitSol(){

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

double DGSolver::initSol_legendre_proj(const int &eID,
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

void DGSolver::ComputeExactSolShift(){

    // Preparing shift information:
    //-------------------------------
    double a=0.;

    wave_length_ = grid_->xf - grid_->x0 ;
    a = simdata_->a_wave_;
    exact_sol_shift = (a * simdata_->t_end_ );

    return;
}

double DGSolver::ExactSol_legendre_proj(const int &eID,
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

void DGSolver::UpdateResid(double **Resid_, double **Qn_){

    register int j;

    // Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    double Ql=0.0,Qr=0.0;

    // fixme: Left and right boundary fluxes :

    j=0.0;

    Ql = evalSolution(&Qn_[grid_->Nelem-1][0], 1.0);
    Qr = evalSolution(&Qn_[j][0], -1.0);

    flux_com[j] = Compute_common_flux(Ql,Qr,simdata_->a_wave_
                                      , simdata_->upwind_param_);

    flux_com[grid_->Nfaces-1] = flux_com[j];

    for(j=1; j<grid_->Nfaces-1; j++){

        Ql = evalSolution(&Qn_[j-1][0], 1.0);
        Qr = evalSolution(&Qn_[j][0], -1.0);

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

void DGSolver::UpdateResidOneCell(const int &cellid, double *q_, double *resid_){

    unsigned int j=cellid;

    int k=0;

    double Mkk=0.0;
    double Lk_p1=0.0, Lk_m1=0.0;
    double f_proj_k=0.0;
    double flux_jp1=0.0;  // f_j+1/2
    double flux_jm1=0.0;  // f_j-1/2
    double fact_=0.0;
    double hjj=0.0;
    double mkk=0.0;

    hjj = grid_->h_j[j];

    fact_ = -2.0/hjj;

    flux_jm1 = flux_com[j];
    flux_jp1 = flux_com[j+1];

    for(k=0; k<Ndof; k++){

        mkk = eval_basis_norm_squared(k);
        Mkk = 1./mkk;
        Lk_m1 = eval_basis_poly(-1,k);
        Lk_p1 = eval_basis_poly( 1,k);
        f_proj_k = eval_localflux_proj(q_,k);

        resid_[k] = fact_ * Mkk* ( flux_jp1 * Lk_p1 - flux_jm1 *Lk_m1 - f_proj_k )  ;

    }

    return;
}

double DGSolver::Compute_common_flux(const double &Ql, const double &Qr
                                      , const double &wave_speed
                                      , const double &upwind_Beta_){

    double f_upw=0.0, f_cent=0.0, f_common_=0.0;
    double aa=wave_speed;
    double BB = upwind_Beta_;

    // f_upw = Rusanov(Ql,Qr,a_wave);

    f_upw = 0.5 * ( aa*(Ql+Qr) - fabs(aa) * (Qr-Ql) );

    f_cent = 0.5 * aa * (Ql+Qr) ;

    f_common_ = (BB * f_upw) + ((1.0-BB) * f_cent );

    return f_common_;
}

void DGSolver::Compute_vertex_sol(){

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

void DGSolver::Compute_exact_vertex_sol(){

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

void DGSolver::Compute_projected_exact_sol(){

    register int j; int k=0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qex_proj[j][k] = ExactSol_legendre_proj(j,k,quad_);

    return;
}

double DGSolver::eval_init_sol(const double& xx){

    if(simdata_->wave_form_==0){

        return sin(2*PI*xx);

    }else if(simdata_->wave_form_==1){

        return exp(-simdata_->Gaussian_exponent_*pow(xx,2));
    }else{
        _notImplemented("Wave form is not implemented");
    }
}

double DGSolver::eval_basis_poly(const double& xi_, const int& basis_k_){

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

double DGSolver::eval_basis_norm_squared(const int &basis_k_){

    // this is norm^2

    return 2./(2*basis_k_+1);
}

double DGSolver::evalSolution(const double* q_, const double& xi_pt_){

    double xx=xi_pt_;
    double Q_=0.0;

    int k;

    for(k=0; k<Ndof; k++){

        Q_ += q_[k] * eval_basis_poly(xx,k);
    }

    return Q_;
}

double DGSolver::eval_localflux_proj(const double *q_, const int &basis_k_){

    int k=basis_k_;

    switch (k) {
    case 0:
        return 0.0;

    case 1:
        return ( 2.0 * simdata_->a_wave_ * q_[k-1] );

    case 2:
        return ( 2.0 * simdata_->a_wave_ * q_[k-1] );

    case 3:
        return ( 2.0 * simdata_->a_wave_ * ( q_[k-3] + q_[k-1]) );

    case 4:
        return ( 2.0 * simdata_->a_wave_ * ( q_[k-3] + q_[k-1]) );

    default:
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"polynomial order of %d ",k-1);
        _notImplemented(ss);
        emptyarray(ss);
        break;

        return 1000.0;
    }


}

double DGSolver::ComputePolyError(){

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

double DGSolver::L1_error_projected_sol(){

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

double DGSolver::L2_error_projected_sol(){

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

double DGSolver::L2_error_nodal_disc_sol(){

    // We have two options: either use only
    // the interface nodes or also use the center node

    double L2_error=0.0;

    return L2_error;
}

double DGSolver::L1_error_average_sol(){

    register int j;

    double L1_error=0.0,error=0.0;

    for(j=0; j<grid_->Nelem; j++){

        error += fabs(Qex_proj[j][0] - Qn[j][0]);
    }

    L1_error = error/grid_->Nelem;

    return L1_error;
}

double DGSolver::L2_error_average_sol(){

    register int j;

    double L2_error=0.0,error=0.0;

    for(j=0; j<grid_->Nelem; j++){

        error += pow((Qex_proj[j][0] - Qn[j][0]),2);
    }

    L2_error = sqrt(error/grid_->Nelem);

    return L2_error;
}

void DGSolver::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[200];

    sprintf(fname,"%snodal/u_nodal_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir, simdata_->Nelem_
            ,simdata_->CFL_,simdata_->upwind_param_
            ,simdata_->Nperiods);

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nfaces; j++)
        fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qv[j]);

    fclose(sol_out);

    emptyarray(fname);

    fname = new char[200];

    sprintf(fname,"%snodal/u_nodal_exact_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir,simdata_->Nelem_
            ,simdata_->CFL_,simdata_->upwind_param_
            ,simdata_->Nperiods);

    FILE* sol_out1=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out1, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact[j]);

    fclose(sol_out1);

    emptyarray(fname);

    return;
}

void DGSolver::print_average_sol(){

    register int j;

    char *fname=nullptr;
    fname = new char[200];

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    sprintf(fname,"%saver/u_aver_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir,simdata_->Nelem_
            ,simdata_->CFL_,simdata_->upwind_param_
            ,simdata_->Nperiods);

    FILE* sol_out=fopen(fname,"w");

     for(j=0; j<grid_->Nelem; j++)
         fprintf(sol_out, "%2.10e %2.10e %2.10e\n"
                 ,grid_->Xc[j], Qn[j][0], Qex_proj[j][0]);

     fclose(sol_out);

     emptyarray(fname);

    return;
}

void DGSolver::dump_errors(double& L1_proj_sol_, double &L1_aver_sol_
                           ,double& L2_proj_sol_, double &L2_aver_sol_){

    char *fname=nullptr;
    fname = new char[200];

    if(simdata_->Sim_mode=="error_analysis_CFL"){
        sprintf(fname,"%serrors/errors_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,simdata_->CFL_
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%d %2.10f %2.10e %2.10e %2.10e %2.10e\n"
                ,grid_->Nelem, time_step
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_);

        fclose(solerror_out);
        emptyarray(fname);

    }else if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%serrors/errors_dt%1.3e_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,time_step,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%d %2.10f %2.10e %2.10e %2.10e %2.10e\n"
                ,grid_->Nelem, CFL
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_);

        fclose(solerror_out);
        emptyarray(fname);

        // Dumping all errors in one file as a function of dt:
        //--------------------------------------------------------
        fname = new char[200];

        sprintf(fname,"%serrors/errors_alldt_N%d_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%1.7e %2.10e %2.10e %2.10e %2.10e\n"
                ,time_step
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_);

         fclose(solerror_out);

         emptyarray(fname);

    }else{
        sprintf(fname,"%serrors/sol_errors_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,simdata_->CFL_,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"w");

        fprintf(solerror_out, "%d %2.10e %2.10e %2.10e %2.10e\n"
                ,grid_->Nelem
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_);

        fclose(solerror_out);

        emptyarray(fname);

        // Dumping all errors in one file as a function of beta:
        //--------------------------------------------------------
        fname = new char[200];

        sprintf(fname,"%serrors/sol_errors_N%d_CFL%1.3f_allBeta_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,simdata_->CFL_
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%1.2f %2.10e %2.10e %2.10e %2.10e\n"
                ,simdata_->upwind_param_
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_);

        fclose(solerror_out);

        emptyarray(fname);
    }

    return;
}

void DGSolver::dump_discont_sol(){

    register int j; int k;

    double xx=0.0,qq=0.0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    char *fname=nullptr;
    fname = new char[200];

    sprintf(fname,"%snodal/u_disc_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,simdata_->Nelem_,simdata_->CFL_
            ,simdata_->upwind_param_,simdata_->Nperiods);

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

    fname = new char[200];

    sprintf(fname,"%snodal/uex_disc_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,simdata_->Nelem_,simdata_->CFL_
            ,simdata_->upwind_param_,simdata_->Nperiods);

    sol_out=fopen(fname,"w");

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


























