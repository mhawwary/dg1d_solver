#include"ExplicitTimeSolver.hpp"

ExplicitTimeSolver::ExplicitTimeSolver(){

 return;
}

ExplicitTimeSolver::~ExplicitTimeSolver(){

    Reset_time_solver();

 return;
}

void ExplicitTimeSolver::setupTimeSolver(DGSolver *dg_solver_, SimData *simdata_){

    space_solver = dg_solver_;
    simdata = simdata_;

    Ndof = space_solver->GetNdof();

    resid = new double* [simdata->Nelem_];

    register int i;

    for(i=0; i<simdata->Nelem_; i++){

        resid[i] = new double[Ndof];
    }

    dt_ = space_solver->GetTimeStep();

    if(simdata_->RK_order_>1){

        q_temp = new double* [simdata->Nelem_];

        for(i=0; i<simdata->Nelem_; i++){

            q_temp[i] = new double[Ndof];
        }
    }

    return;
}

void ExplicitTimeSolver::SolveOneStep(double **qn_){

    space_solver->UpdateResid(resid,qn_);

    switch (simdata->RK_order_) {

    case 1:

        FwdEuler(qn_);

        break;

    case 2:

        SSPRK22(qn_);

        break;

    case 3:

        SSPRK33(qn_);

        break;

    default:
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"RK order of %d ",simdata->RK_order_);
        _notImplemented(ss);
        emptyarray(ss);
        break;
    }

    IterNo++;

    return;
}

void ExplicitTimeSolver::CopyOldSol(double **q_t_, double **qn_){

    register int j;

    unsigned int k;

    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j][k];

    return;
}

void ExplicitTimeSolver::FwdEuler(double **qn_){

    register int i;
    unsigned int j;


    for(i=0; i<simdata->Nelem_; i++)
        for(j=0; j<Ndof; j++)
            qn_[i][j] = qn_[i][j] + dt_ * resid[i][j];

    return;
}

void ExplicitTimeSolver::SSPRK22(double **q_){

    register int j;

    unsigned int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it

    // Step1:
    //-----------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++){
            q_[j][k] = q_temp[j][k] + dt_ * resid[j][k];
        }

    space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++) {
            q_[j][k] = 0.5 * ( q_temp[j][k] +  q_[j][k]
                                      + dt_ * resid[j][k] );
        }

    return;
}

void ExplicitTimeSolver::SSPRK33(double **q_){

    register int j;

    unsigned int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it


    // Step1:
    //-----------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] + dt_ * resid[j][k];


    space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] =  (0.75 * q_temp[j][k] )
                    + 0.25 * ( q_[j][k] + dt_ * resid[j][k] );

    space_solver->UpdateResid(resid,q_);

    // Step3:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] =  ( q_temp[j][k]/3. )
                    + 2. * ( q_[j][k] + dt_ * resid[j][k] ) / 3.;

    return;
}


void ExplicitTimeSolver::Reset_time_solver(){

    emptyarray(simdata->Nelem_,resid);
    emptyarray(simdata->Nelem_,q_temp);

    //emptypointer(space_solver);
    //emptypointer(simdata);

    _print("finshed deallocating time solver");

    return;
}


void ExplicitTimeSolver::UpdateIter(){

    IterNo++;

    return;
}

int ExplicitTimeSolver::GetIter(){

    return IterNo;
}
