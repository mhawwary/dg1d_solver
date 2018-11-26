#include "DGSolverAdvec.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

DGSolverAdvec::~DGSolverAdvec(void){

    Reset_solver();
}

void DGSolverAdvec::setup_solver(GridData& meshdata_, SimData& osimdata_){

    register int i;

    grid_ = &meshdata_;
    simdata_ = &osimdata_;

    Ndof= simdata_->poly_order_+1;

    // Nquad  is for error integrations
    // Nquad_invFlux is for flux projection
    switch (Ndof) {
    case 1:  // p0
        Nquad_ = 8;  // 1, but use highest possible for more accuracy
        Nquad_invFlux_= 1; // it is not applicable
        break;
    case 2: // p1
        Nquad_ = 8;  // 2, but use highest possible for more accuracy
        Nquad_invFlux_= 2;
        break;
    case 3:  // p2
        Nquad_ = 8;  // 3, but use highest possible for more accuracy
        Nquad_invFlux_= 3;
        break;
    case 4:  // p3
        Nquad_ = 8;  // 4, but use highest possible for more accuracy
        Nquad_invFlux_= 5;
        break;
    case 5:  // p4
        Nquad_ = 8;  // 5, but use highest possible for more accuracy
        Nquad_invFlux_= 6;
        break;
    case 6:  // p5
        Nquad_ = 8;  // 6, but use highest possible for more accuracy
        Nquad_invFlux_= 8;
        break;
    default:
        break;
    } // This needs to be tested in the future

    quad_.setup_quadrature(Nquad_);
    quad_invF_.setup_quadrature(Nquad_invFlux_);

    Lk = new double*[Ndof]; // Lk[k, xi]
    Lk_norm_squar = new double[Ndof];

    for(i=0; i<Ndof; i++)
        Lk[i] = new double[2];

    Qn    =  new double* [grid_->Nelem];
    Qex_proj = new double*[grid_->Nelem];

    for(i=0; i<grid_->Nelem; i++){
        Qn[i]       = new double[Ndof];
        Qex_proj[i] = new double[Ndof];
    }

    Q_exact = new double[grid_->N_exact_ppts];
    Qv = new double[grid_->Nfaces];

    flux_com = new double[grid_->Nfaces];
    Q_cont_sol = new double[grid_->N_uniform_pts];

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

void DGSolverAdvec::setup_basis_interpolation_matrices(){
    int k;
    for(k=0; k<Ndof; k++){
        Lk[k][0] = eval_basis_poly(-1,k);
        Lk[k][1] = eval_basis_poly(1,k);
        Lk_norm_squar[k] = eval_basis_norm_squared(k);
    }
    return;
}

void DGSolverAdvec::Reset_solver(){

    emptyarray(grid_->Nelem,Qn);
    emptyarray(Q_exact);
    //emptyarray(qq_exact_time);
    emptyarray(flux_com);
    emptyarray(Qv);
    emptyarray(grid_->Nelem,Qex_proj);

    emptyarray(Q_cont_sol);

    emptyarray(Ndof,Lk);
    emptyarray(Lk_norm_squar);

    quad_.Reset_quad();
    quad_invF_.Reset_quad();

    grid_->Reset_();

    printf("\nMay need to add simdata reset like the other solvers ???\n");
    // May need to add simdata reset like the other solvers ???
    return;
}


// Solver functions
//-------------------------------------------

void DGSolverAdvec::CalcTimeStep(){

    compute_uniform_cont_sol();
    double TV_ = compute_totalVariation();

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;

    if(simdata_->calc_dt_flag==1){  // use CFL as input
        CFL = simdata_->CFL_;
        time_step = (grid_->dx * CFL )/ simdata_->a_wave_;
        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){  // use dt as input
        time_step = simdata_->dt_;
        last_time_step = time_step;
        CFL = simdata_->a_wave_ * time_step / grid_->dx ;
        simdata_->CFL_ = CFL;

    }else {

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
    cout << "CFL no.        : "<<CFL<<endl;
    cout << "time step, dt  : "<<time_step<<endl;
    cout << "last_time_step : "<<last_time_step<<endl;
    cout << "input Nperiods : "<<simdata_->Nperiods<<endl;
    cout << "new   Nperiods : "<<simdata_->t_end_/T_period<<endl;
    cout << "exact_sol_shift: "<<wave_speed_*simdata_->a_wave_<<endl;
    cout << "T_period       : "<<T_period<<endl;
    printf("actual_end_time:%1.5f",simdata_->t_end_);
    cout <<"\nMax_iter: "<<simdata_->maxIter_<<endl;

    cout << "\nNumber of Elements: "<< grid_->Nelem<<"  dx:  "<<grid_->dx<<endl;
    cout << "Polynomial  order : "<< simdata_->poly_order_  << endl;
    cout << "Runge-Kutta order : "<< simdata_->RK_order_    << endl;
    cout << "Upwind parameter  : "<< simdata_->upwind_param_<< endl;
    cout << "Poly GaussQuad order  : "<< Nquad_ << endl;
    cout << "Flux GaussQuad order  : "<< Nquad_invFlux_ << endl;
    cout <<"===============================================\n";

    return;
}

void DGSolverAdvec::InitSol(){

    register int j;
    int k=0,i;
    double xx=0.0,qi_=0.0;

    max_eigen_advec=0.0;  // initializing the maximum eigen value with zero
    GaussQuad quad_temp_; quad_temp_.setup_quadrature(8);

    if(simdata_->eqn_type_=="inv_burger"){  // burger's equation
        max_eigen_advec=0.0;  // initializing the maximum eigen value with zero
        for(j=0; j<grid_->Nelem; j++){
            for(k=0; k<Ndof; k++)
                Qn[j][k] = initSol_legendre_proj(j,k,quad_temp_);

            for (i=0; i<quad_temp_.Nq; i++){
                xx = 0.5 * grid_->h_j[j] * quad_temp_.Gaus_pts[i] + grid_->Xc[j];
                qi_ = evalSolution(&Qn[j][0],xx);
                if(fabs(qi_)>max_eigen_advec) max_eigen_advec = fabs(qi_);
            }
        }

    }else if(simdata_->eqn_type_=="linear_advec"){ // linear wave equation
        for(j=0; j<grid_->Nelem; j++)
            for(k=0; k<Ndof; k++)
                Qn[j][k] = initSol_legendre_proj(j,k,quad_temp_);

        max_eigen_advec = simdata_->a_wave_;

    }else{
        FatalError_exit("Equation type is not implemented");
    }

    CalcTimeStep(); // based on maximum eigenvalues
    quad_temp_.Reset_quad();

    init_wave_E_ = Compute_waveEnergy(Qn);  // initial wave energy of the projected solution

    return;
}

double DGSolverAdvec::initSol_legendre_proj(const int &eID,
                                            const int &basis_k,
                                            const GaussQuad &quad_){
    int i=0,j=0,k=0;
    k=basis_k;
    j=eID;
    double xx=0.0;
    double II=0.0;
    double Qinit_=0.0;
    double Lk_=1.0;

    II=0.0;
    for (i=0; i<quad_.Nq; i++){
        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
        Qinit_= eval_init_sol(xx);
        Lk_ = eval_basis_poly(quad_.Gaus_pts[i], k);
        II += quad_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverAdvec::ComputeExactSolShift(){
    exact_sol_shift = (wave_speed_ * phy_time );
    return;
}

double DGSolverAdvec::ExactSol_legendre_proj(const int &eID,
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

        if(simdata_->wave_form_==0){    // single mode wave
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

double DGSolverAdvec::TimeAccurateExactSol_legendre_proj(const int &eID,
                                                         const int &basis_k,
                                                         const GaussQuad &quad_temp_){
    int i=0,j=0,k=0;
    k=basis_k;
    j=eID;
    double xx=0.0;
    double x0,x1;
    double II=0.0;
    double Qinit_=0.0;
    double Lk_=1.0;
    double time_accurate_shift=0.0;
    double a_wave_=0.0;
    a_wave_ = simdata_->a_wave_;
    time_accurate_shift =  (a_wave_ * phy_time);

    //printf("Physical Time= %1.5f\n",phy_time);
    //std::cin.get();

    II=0.0;
    for (i=0; i<quad_temp_.Nq; i++){
        xx = 0.5 * grid_->h_j[j] * quad_temp_.Gaus_pts[i]
                + grid_->Xc[j] - time_accurate_shift;

        if(simdata_->wave_form_==0){    // single mode wave
            Qinit_ = eval_init_sol(xx);

        }else if(simdata_->wave_form_==1){  // Gaussian wave
            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);
            Qinit_= eval_init_sol(x0)+eval_init_sol(x1);
            if(x0==0 && x1==0) Qinit_ = 0.5*Qinit_;
        }
        Lk_ = eval_basis_poly(quad_temp_.Gaus_pts[i], k);
        II += quad_temp_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }
    II = II / Lk_norm_squar[k] ;

    return II;
}

void DGSolverAdvec::UpdateResid(double **Resid_, double **Qn_){

    register int j;
    double Ql=0.0,Qr=0.0;

    // Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    // fixme: Left and right boundary fluxes :
    j=0.0;
    Ql = evalSolution(&Qn_[grid_->Nelem-1][0], 1);
    Qr = evalSolution(&Qn_[j][0], 0);
    flux_com[j] = Compute_common_flux(Ql,Qr,simdata_->a_wave_
                                      , simdata_->upwind_param_);
    flux_com[grid_->Nfaces-1] = flux_com[j];

    // Interior Faces:
    for(j=1; j<grid_->Nfaces-1; j++){
        Ql = evalSolution(&Qn_[j-1][0], 1);
        Qr = evalSolution(&Qn_[j][0], 0);
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

void DGSolverAdvec::UpdateResidOneCell(const int &cellid, double *q_, double *resid_){

    // General parameters:
    unsigned int j=cellid;
    int k=0;
    double Mkk=0.0;
    double Lk_p1=0.0, Lk_m1=0.0;
    double fact_=0.0;
    double hjj=0.0;
    hjj = grid_->h_j[j];
    fact_ = 2.0/hjj;

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
        f_proj_k = eval_local_invflux_proj_exact(q_,k);
        resid_[k] = fact_ * Mkk
                * (- ( flux_jp1 * Lk_p1 - flux_jm1 *Lk_m1) + f_proj_k );
    }

    return;
}

double DGSolverAdvec::Compute_common_flux(const double &Ql, const double &Qr
                                          , const double &wave_speed
                                          , const double &upwind_Beta_){

    double f_upw=0.0, f_cent=0.0, f_common_=0.0;
    double aa=wave_speed;
    double BB = upwind_Beta_;
    double Fl=0.0,Fr=0.0;

    if(simdata_->eqn_type_=="inv_burger"){ // burger's equation
        f_upw = Rusanov(Ql,Qr);
        Fl = 0.5 * pow(Ql,2);
        Fr = 0.5 * pow(Qr,2);

    }else if(simdata_->eqn_type_=="linear_advec"){  // linear wave equation
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

double DGSolverAdvec::Rusanov(const double &Ql, const double &Qr){

    // Now it is only working for burger's equation:
    double Lambda_max = 0.0;
    double Fl=0.0,Fr=0.0;

    Fl = 0.5 *  pow(Ql,2);
    Fr = 0.5 *  pow(Qr,2);
    Lambda_max = 0.5 * fabs(Ql+Qr);

    return  ( 0.5 * (Fl+Fr) - 0.5 * Lambda_max * (Qr-Ql) );
}

double DGSolverAdvec::eval_burgers_invflux(const double& xi_pt_
                                           , const double *q_){
    int k;
    double xx=xi_pt_;
    double Q_;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_basis_poly(xx,k);

    return ( 0.5 * pow(Q_,2) );
}

double DGSolverAdvec::eval_local_invflux_proj_exact(const double *q_, const int &basis_k_){

    int i;
    double II=0.0;

    if(simdata_->eqn_type_=="inv_burger"){   // burger's equation
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

    }else if(simdata_->eqn_type_=="linear_advec"){    // linear wave equation
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

void DGSolverAdvec::Compute_vertex_sol(){

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

void DGSolverAdvec::Compute_exact_vertex_sol(){

    register int j;
    double xx=0.0;
    double x0,x1;

    ComputeExactSolShift(); // updating the shift for the currnet time

    for(j=0; j<grid_->N_exact_ppts; j++){
        xx = grid_->x_exact_ppts[j]- exact_sol_shift;

        if(simdata_->wave_form_==0){   // single mode wave
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

void DGSolverAdvec::Compute_TimeAccurate_exact_sol(){

    register int j;
    double xx=0.0;
    double x0,x1;
    double time_accurate_shift=0.0;
    double a_wave_=0.0;
    a_wave_ = simdata_->a_wave_;
    time_accurate_shift = a_wave_*phy_time;

    for(j=0; j<grid_->N_exact_ppts; j++){
        xx = grid_->x_exact_ppts[j]- time_accurate_shift;

        if(simdata_->wave_form_==0){   // single mode wave
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

void DGSolverAdvec::Compute_projected_exact_sol(){
    register int j; int k=0;
    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qex_proj[j][k] = ExactSol_legendre_proj(j,k,quad_);

    return;
}

double DGSolverAdvec::eval_init_sol(const double& xx){

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

double DGSolverAdvec::eval_exact_sol(double &xx){

    //xx = xx - exact_sol_shift;
    double time_accurate_shift = simdata_->a_wave_ * phy_time;
    xx = xx - time_accurate_shift;

    if(simdata_->wave_form_==0){  // single mode wave
        return eval_init_sol(xx);

    }else if(simdata_->wave_form_==1){ // Gaussian wave
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

double DGSolverAdvec::eval_basis_poly(const double& xi_, const int& basis_k_){

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

    }else if(basis_k_==5) {
        return  (63.0 * pow(xi_,5) - 70.0 * pow(xi_,3) + 15.0 * xi_ ) /8. ;

    }else {
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"polynomial order of: %d",basis_k_);
        _notImplemented(ss);

        return 0.0;
    }
}

double DGSolverAdvec::eval_basis_norm_squared(const int &basis_k_){
    // this is norm^2
    return 2./(2*basis_k_+1);
}

double DGSolverAdvec::evalSolution(const double* q_, const double& xi_pt_){

    double xx=xi_pt_;
    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_basis_poly(xx,k);

    return Q_;
}

double DGSolverAdvec::evalSolution(const double* q_, const int& i_pos_){

    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * Lk[k][i_pos_];   // i_pos_ = 0 if xi =-1, i_pos_ =1 if xi =1

    return Q_;
}

double DGSolverAdvec::ComputePolyError(){

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

double DGSolverAdvec::L1_error_nodal_cont_sol(){

    register int j;
    double L1_error=0.0,II=0.0;
    II=0.0;

    for(j=0; j<grid_->N_uniform_pts; j++)
        II += fabs(Q_exact[j] - Q_cont_sol[j]);

    L1_error = II/grid_->N_uniform_pts;

    return L1_error;
}

double DGSolverAdvec::L2_error_nodal_cont_sol(){

    register int j;
    double L2_error=0.0,II=0.0;
    II=0.0;

    for(j=0; j<grid_->N_uniform_pts; j++)
        II += pow((Q_exact[j] - Q_cont_sol[j]),2);

    L2_error = sqrt(II/grid_->N_uniform_pts);

    return L2_error;
}

double DGSolverAdvec::L1_error_nodal_gausspts_proj(){

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

double DGSolverAdvec::L2_error_nodal_gausspts_proj(){

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

double DGSolverAdvec::L1_error_nodal_gausspts(){

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

double DGSolverAdvec::L2_error_nodal_gausspts(){

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

double DGSolverAdvec::L1_error_projected_sol(){

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

double DGSolverAdvec::L2_error_projected_sol(){

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

double DGSolverAdvec::L1_error_average_sol(){

    register int j;
    double L1_error=0.0,error=0.0;
    for(j=0; j<grid_->Nelem; j++)
        error += fabs(Qex_proj[j][0] - Qn[j][0]);
    L1_error = error/grid_->Nelem;

    return L1_error;
}

double DGSolverAdvec::L2_error_average_sol(){
    register int j;
    double L2_error=0.0,error=0.0;
    for(j=0; j<grid_->Nelem; j++)
        error += pow((Qex_proj[j][0] - Qn[j][0]),2);
    L2_error = sqrt(error/grid_->Nelem);

    return L2_error;
}

void DGSolverAdvec::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[100];

    if(simdata_->Sim_mode=="error_analysis_dt"){
        sprintf(fname,"%snodal/u_cont_N%d_dt%1.3e_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nfaces; j++)
            fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qv[j]);

        fclose(sol_out);
        emptyarray(fname);

    } else{

        sprintf(fname,"%snodal/u_cont_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
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

void DGSolverAdvec::print_average_sol(){

    register int j;

    char *fname=nullptr;
    fname = new char[100];

    if(simdata_->Sim_mode=="error_analysis_dt"){

        sprintf(fname,"%saver/u_aver_N%d_dt%1.3e_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++)
            fprintf(sol_out, "%2.10e %2.10e %2.10e\n"
                    ,grid_->Xc[j], Qex_proj[j][0], Qn[j][0]);

        fclose(sol_out);
        emptyarray(fname);

    }else{

        sprintf(fname,"%saver/u_aver_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
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

void DGSolverAdvec::dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
                                ,double& L1_aver_sol_,double& L2_aver_sol_
                                ,double& L1_nodal_gausspts, double& L2_nodal_gausspts){
    char *fname=nullptr;
    fname = new char[100];

    if(simdata_->Sim_mode=="error_analysis_CFL"){
        sprintf(fname,"%serrors/errors_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,CFL
                ,simdata_->upwind_param_
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

        sprintf(fname,"%serrors/errors_N%d_dt%1.3e_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
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

        sprintf(fname,"%serrors/errors_N%d_alldt_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,simdata_->upwind_param_
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

        sprintf(fname,"%serrors/errors_dt%1.3e_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,time_step
                ,simdata_->upwind_param_
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
        sprintf(fname,"%serrors/errors_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"w");

        fprintf(solerror_out, "%2.10e %2.10e %2.10e %2.10e %2.10e %2.10e\n"
                ,L1_proj_sol_, L1_aver_sol_
                ,L2_proj_sol_, L2_aver_sol_
                ,L1_nodal_gausspts,L2_nodal_gausspts);

        fclose(solerror_out);

        emptyarray(fname);

    }else if(simdata_->Sim_mode=="error_analysis_Beta"){

        sprintf(fname,"%serrors/errors_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
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

        sprintf(fname,"%serrors/errors_N%d_CFL%1.3f_allBeta_%1.3fT.dat"
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

void DGSolverAdvec::dump_timeaccurate_errors(){

    //Fix me! You always need to make sure that the exact solution has been updated, i.e,
    // Compute_projected_exact_sol(); Compute_exact_vertex_sol(); has been called
    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"%serrors/errors_N%d_CFL%1.4f_Beta%1.2f_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,CFL
            ,simdata_->upwind_param_
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

void DGSolverAdvec::dump_discont_sol(){

    register int j; int k;
    double xx=0.0,qq=0.0;
    char *fname=nullptr;
    fname = new char[100];

    if(simdata_->Sim_mode=="error_analysis_dt"){
        sprintf(fname,"%snodal/u_disc_N%d_dt%1.3e_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++)
            for(k=0; k<grid_->N_xi_disc_ppts; k++) {
                qq = evalSolution(&Qn[j][0],grid_->xi_disc[k]);
                xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                        + grid_->Xc[j];
                fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
            }

        fclose(sol_out);
        emptyarray(fname);

    }else {
        sprintf(fname,"%snodal/u_disc_N%d_CFL%1.3f_Beta%1.2f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<grid_->Nelem; j++)
            for(k=0; k<grid_->N_xi_disc_ppts; k++) {
                qq = evalSolution(&Qn[j][0],grid_->xi_disc[k]);
                xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                        + grid_->Xc[j];
                fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
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

    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<grid_->N_xi_disc_ppts; k++) {
            qq = evalSolution(&Qex_proj[j][0],grid_->xi_disc[k]);
            xx = ( 0.5 * grid_->h_j[j] * grid_->xi_disc[k])
                    + grid_->Xc[j];

            fprintf(sol_out,"%2.10e %2.10e\n",xx,qq);
        }

    fclose(sol_out);
    emptyarray(fname);

    return;
}

void DGSolverAdvec::dump_timeaccurate_sol(){

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
        sprintf(fname,"%stime_data/u_cont_N%d_CFL%1.4f_Beta%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        // Dump time accurate continuous equally spaced solution data:
        sprintf(fname,"%stime_data/u_cont_N%d_dt%1.3e_Beta%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
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
        sprintf(fname,"%stime_data/u_disc_N%d_CFL%1.4f_Beta%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->upwind_param_
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        sprintf(fname,"%stime_data/u_disc_N%d_dt%1.3e_Beta%1.2f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->upwind_param_
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

    //Continuous Exact solution:
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

void DGSolverAdvec::compute_uniform_cont_sol(){

    double qq=0.0,qq_end_=0.0;
    register int j; int k;

    // Dump continuous data on uniform points:
    //-------------------------------------------
    int count_=0; //continuous points counter

    // For element zero:
    k = simdata_->N_uniform_pts_per_elem_-1;
    j= grid_->Nelem-1;
    qq = evalSolution(&Qn[j][0],grid_->xi_uniform[k]); // eval last element solution

    k=0; j=0;
    qq += evalSolution(&Qn[j][0],grid_->xi_uniform[k]); // eval first element solution
    qq = 0.5*qq; qq_end_=qq;
    Q_cont_sol[count_] = qq;
    count_++;

    for(k=1; k<simdata_->N_uniform_pts_per_elem_-1; k++) {
        qq = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
        Q_cont_sol[count_] = qq;
        count_++;
    }

    qq=0.0;
    for(j=1; j<grid_->Nelem; j++){

        k=0;
        qq = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
        k = simdata_->N_uniform_pts_per_elem_-1;
        qq += evalSolution(&Qn[j-1][0],grid_->xi_uniform[k]);

        qq = 0.5*qq;
        Q_cont_sol[count_] = qq;
        count_++;
        qq=0.0;

        for(k=1; k<simdata_->N_uniform_pts_per_elem_-1; k++) {
            qq = evalSolution(&Qn[j][0],grid_->xi_uniform[k]);
            Q_cont_sol[count_] = qq;
            count_++;
        }
    }
    Q_cont_sol[count_] = qq_end_;

    //printf("\n Count: %d \t N_uniform: %d\n", count_, grid_->N_uniform_pts);

    return;
}

double DGSolverAdvec::compute_totalVariation(){

    register int i;
    double TV_=0.0;
    for(i=1; i<grid_->N_uniform_pts; i++)
        TV_ += fabs(Q_cont_sol[i]-Q_cont_sol[i-1]);

    return TV_;
}




















