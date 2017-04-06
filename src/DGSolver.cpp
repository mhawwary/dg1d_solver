#include "DGSolver.hpp"


DGSolver::DGSolver(void){

    //Reset_solver();

}

DGSolver::~DGSolver(void){

    Reset_solver();
}

void DGSolver::setup_solver(GridData *meshdata_, SimData* osimdata_){

    grid_ = new GridData;
    simdata_ = new SimData;

    //grid_->set_grid_param(simdata_);

    cout << "\nfinished setuping simdata_ and grid_ in DGSolver\n";

    grid_ = meshdata_;
    simdata_ = osimdata_;

    cout << "\nfinished copying simdata_ and grid_ in DGSolver\n";

    //Nelem = grid_->Nelem;

    //Nfaces = grid_->Nfaces;

    //X = grid_->X;

    //h_j = grid_->h_j;

    //poly_order = simdata_->poly_order_;

    Ndof= simdata_->poly_order_+1;

    Qn    =  new double* [grid_->Nelem];
    Resid =  new double* [grid_->Nelem];

    register int i;

    for(i=0; i<grid_->Nelem; i++){

        Qn[i]    = new double[Ndof];
        Resid[i] = new double[Ndof];
    }

    Q_exact = new double[grid_->Nfaces];

    flux_com = new double[grid_->Nfaces];

    //Beta_ = simdata_.upwind_param_;



    return;
}

void DGSolver::InitSol(){

    register int i,k,j;

    GaussQuad quad_;

    quad_.setup_quadrature(5);

    for(j=0; j<grid_->Nelem; j++){

        for(k=0; k<Ndof; k++){

            Qn[j][k] = initSol_legendre_proj(j,k,quad_);
        }
    }

    return;
}

double DGSolver::eval_init_sol(const double& xx){

    return sin(2*PI*xx);
}

double DGSolver::initSol_legendre_proj(const int &eID,
                                       const int &basis_id,
                                       const GaussQuad &quad_){

    register int i; int k=0,j=0;

    double xx=0.0;
    double II=0.0;
    double Qinit_=0.0;
    double _phy_=1.0;
    double _phy_norm=1.0;

    k=basis_id;
    j=eID;

    II=0.0;

    for (i=0; i<quad_.Nq; i++){

        xx = 0.5 * grid_->h_j[j] * quad_.Gaus_pts[i] + grid_->Xc[j];
        Qinit_= eval_init_sol(xx);

        //_phy_ = eval_basis_poly(quad_.Gaus_pts[i], k);

        II += Qinit_ * _phy_ ;
    }

    II = II/(_phy_norm*_phy_norm);

    return II;
}

void DGSolver::Reset_solver(){

    //emptyarray(X);
    //emptyarray(Xc);
    //emptyarray(h_j);
    emptyarray(xi);

    emptyarray(grid_->Nelem,Qn);
    emptyarray(Q_exact);
    emptyarray(flux_com);
    emptyarray(grid_->Nelem,Resid);

    emptyarray(Mnn);
    emptyarray(phy_n);

    return;
}


















