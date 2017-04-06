#include "DGSolver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
DGSolver::DGSolver(void){

    //Reset_solver();

}

DGSolver::~DGSolver(void){

    Reset_solver();
}

void DGSolver::setup_solver(GridData& meshdata_, SimData& osimdata_){

    grid_ = new GridData;
    simdata_ = new SimData;

    //grid_->set_grid_param(simdata_);

    cout << "\n--finished setuping simdata_ and grid_ in DGSolver\n";

    grid_ = &meshdata_;
    simdata_ = &osimdata_;

    cout << "\n--finished copying simdata_ and grid_ in DGSolver\n";


    Ndof= simdata_->poly_order_+1;

    Qn    =  new double* [grid_->Nelem];
    //Resid =  new double* [grid_->Nelem];

    register int i;

    for(i=0; i<grid_->Nelem; i++){

        Qn[i]    = new double[Ndof];
        //Resid[i] = new double[Ndof];
    }

    Q_exact = new double[grid_->Nfaces];
    Qv = new double[grid_->Nfaces];

    flux_com = new double[grid_->Nfaces];

    SetPhyTime(simdata_->t_init_);

    CalcTimeStep();

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

    grid_->Reset_();

//    _print("finshed Reseting grid_ in DGsolver");

//    emptypointer(grid_);

//    _print("finshed deallocating grid_ in DG solver");

//    emptypointer(simdata_);

//    _print("finshed deallocating DG solver");


    return;
}


// Solver functions
//-------------------------------------------

void DGSolver::SetPhyTime(const double &time_){

    phy_time=time_;

    return;
}

void DGSolver::CalcTimeStep(){

    if(simdata_->calc_dt_flag==1){

        time_step = (grid_->dx * simdata_->CFL_ )/ simdata_->a_wave_;


    }else if(simdata_->calc_dt_flag==0){

        time_step = simdata_->dt_;

        CFL = simdata_->a_wave_ * time_step / grid_->dx ;

    }else {

        FatalError("Wrong Calc_dt_flag");
    }

    return;
}

double DGSolver::GetTimeStep(){

    return time_step;
}

double DGSolver::GetCFL(){

    return CFL;
}

unsigned int DGSolver::GetNdof(){

    return Ndof;
}

double** DGSolver::GetNumSolution(){

    return Qn;
}

double* DGSolver::GetVertexNumSol(){

    return Qv;
}

double* DGSolver::GetExactSolution(){

    return Q_exact;
}

double DGSolver::GetPhyTime(){

    return phy_time;
}


void DGSolver::InitSol(){

    register int j;

    unsigned int k=0;


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

    unsigned int i;
    unsigned int k=0;

    int j=0;

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

    /*if(k==0){

       double qq=(cos(2*PI* grid_->X[j])-cos(2*PI*grid_->X[j+1] ))/(2*PI*grid_->h_j[j]);
        _compare(II,qq);
        if(fabs(II-qq)>=1e-8)   FatalError("\nInitial sol Gauss projection is wrong");

    }else if(k==1){

        double qq= (3./(2*PI*PI*pow(grid_->h_j[j],2))) * ( sin(2*PI*grid_->X[j+1])-sin(2*PI*grid_->X[j] )
                - PI*grid_->h_j[j]* ( cos(2*PI*grid_->X[j+1]) + cos(2*PI*grid_->X[j]) ) );
         _compare(II,qq);
         if(fabs(II-qq)>=1e-8)   FatalError("\nInitial sol Gauss projection is wrong");
    }*/

    return II;
}

void DGSolver::UpdatePhyTime(const double& dt_){

    phy_time += dt_;

    return;
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

    //_print("finished computing common fluxes inside space solver");

    // Element loop to calculate and update the residual:
    //----------------------------------------------------

    for(j=0; j<grid_->Nelem; j++){
         UpdateResidOneCell(j, &Qn_[j][0], &Resid_[j][0]);
    }

    return;
}

void DGSolver::UpdateResidOneCell(const int &cellid, double *q_, double *resid_){

    unsigned int j=cellid;

    unsigned int k=0;

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

    f_upw = 0.5 * (aa*(Ql+Qr) - fabs(aa) * (Qr-Ql) );

    f_cent = 0.5 * aa * (Ql+Qr) ;

    f_common_ = (BB * f_upw + (1.0-BB) * f_cent );

    //_print("finished computing common flux inside commonflux function");

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

void DGSolver::Compute_exact_sol(){

    register int j;

    for(j=0; j<grid_->Nfaces; j++){

        Q_exact[j] = eval_init_sol(grid_->X[j]-simdata_->a_wave_*simdata_->t_end_);
    }

    return;
}

double DGSolver::eval_init_sol(const double& xx){

    return sin(2*PI*xx);
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

    unsigned int k;

    for(k=0; k<Ndof; k++){

        Q_ += q_[k] * eval_basis_poly(xx,k);
    }

    return Q_;
}

double DGSolver::eval_localflux_proj(const double *q_, const unsigned int &basis_k_){

    unsigned int k=basis_k_;

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


void DGSolver::print_num_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./output/u_final.dat");

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nfaces; j++){

        fprintf(sol_out, "%2.10e\n", Qv[j]);
    }

    fclose(sol_out);

    emptyarray(fname);

    return;
}


void DGSolver::print_exact_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./output/u_exact.dat");

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nfaces; j++){

        fprintf(sol_out, "%2.10e\n", Q_exact[j]);
    }

    fclose(sol_out);

    emptyarray(fname);


    return;
}


void DGSolver::print_exact_average_sol(){

    register int j;

    double Q_exact_aver_=0.0;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./output/u_aver_exact.dat");

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->Nelem; j++){

        Q_exact_aver_ = initSol_legendre_proj(j,0,quad_);

        fprintf(sol_out, "%2.10e\n", Q_exact_aver_);
    }

    fclose(sol_out);

    emptyarray(fname);

    return;
}


void DGSolver::print_num_average_sol(){

    register int j;

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./output/u_aver_final.dat");

    FILE* sol_out=fopen(fname,"w");

     for(j=0; j<grid_->Nelem; j++){

         fprintf(sol_out, "%2.10e\n", Qn[j][0]);

     }

     fclose(sol_out);

     emptyarray(fname);

     char *fname1=nullptr;
     fname1 = new char[100];

     sprintf(fname1,"./output/Xc.dat");

     FILE* sol_out1=fopen(fname1,"w");

     for(j=0; j<grid_->Nelem; j++){

         fprintf(sol_out1, "%2.10e\n", grid_->Xc[j]);

     }

     fclose(sol_out1);

     emptyarray(fname1);

    return;
}































