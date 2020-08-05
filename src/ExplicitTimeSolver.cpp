#include"ExplicitTimeSolver.hpp"

ExplicitTimeSolver::ExplicitTimeSolver(){

 return;
}

ExplicitTimeSolver::~ExplicitTimeSolver(){

    Reset_time_solver();

 return;
}

void ExplicitTimeSolver::setupTimeSolver(DGSolver *dg_solver_
                                         , SimData *osimdata_){

    space_solver = dg_solver_;

    simdata = osimdata_;

    Ndof = space_solver->GetNdof();

    resid = new double* [simdata->Nelem_];

    register int i;

    for(i=0; i<simdata->Nelem_; i++){

        resid[i] = new double[Ndof];
    }

    dt_ = space_solver->GetTimeStep();

    if(simdata->RK_order_>1){

        q_temp = new double* [simdata->Nelem_];

        for(i=0; i<simdata->Nelem_; i++)
            q_temp[i] = new double[Ndof];

        if(simdata->RK_order_==4){

            resid_temp0 = new double* [simdata->Nelem_];
            resid_temp1 = new double* [simdata->Nelem_];
            resid_temp2 = new double* [simdata->Nelem_];

            for(i=0; i<simdata->Nelem_; i++){
              resid_temp0[i] = new double[Ndof];
              resid_temp1[i] = new double[Ndof];
              resid_temp2[i] = new double[Ndof];
            }
        }else if(simdata->RK_order_==9){
          q_temp1 = new double* [simdata->Nelem_];
          q_temp2 = new double* [simdata->Nelem_];
          q_temp3 = new double* [simdata->Nelem_];
          q_temp4 = new double* [simdata->Nelem_];
          q_temp5 = new double* [simdata->Nelem_];
          q_temp6 = new double* [simdata->Nelem_];
          q_temp7 = new double* [simdata->Nelem_];

          for(i=0; i<simdata->Nelem_; i++){
              q_temp1[i] = new double[Ndof];
              q_temp2[i] = new double[Ndof];
              q_temp3[i] = new double[Ndof];
              q_temp4[i] = new double[Ndof];
              q_temp5[i] = new double[Ndof];
              q_temp6[i] = new double[Ndof];
              q_temp7[i] = new double[Ndof];
          }
        }
    }

    return;
}

void ExplicitTimeSolver::SolveOneStep(double **qn_){

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

    case 4:

        classicRK4(qn_);

        break;

    case 9:

        SSPRK99(qn_);

        break;

    default:
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"RK order of %d ",simdata->RK_order_);
        _notImplemented(ss);
        emptyarray(ss);
        break;
    }

    IterNo++;

//    if(IterNo==simdata->maxIter_-1)
//        dt_ = space_solver->GetLastTimeStep();

    return;
}

void ExplicitTimeSolver::CopyOldSol(double **q_t_, double **qn_){

    register int j;

    int k;

    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j][k];

    return;
}

void ExplicitTimeSolver::CopyOldResid(double **resid_t_, double **old_resid_){

    register int j;

    int k;

    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            resid_t_[j][k] = old_resid_[j][k];

    return;
}

void ExplicitTimeSolver::FwdEuler(double **q_){

    register int i;
    int j;

    for(i=0; i<simdata->Nelem_; i++)
        for(j=0; j<Ndof; j++)
            q_[i][j] = q_[i][j] + dt_ * resid[i][j];

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK22(double **q_){

    register int j;

    int k;

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

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK33(double **q_){

    register int j;

    int k;

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

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK99(double **q_){

    register int j,i;

    int k;

    double coeff[9]={1.0,2119./5760.,103./560.,53./864.,11./720.,1./320.,1./2160.,1./10080.,1./362880.};
    for(i=1;i<9; i++)
      coeff[0]-=coeff[i];

    CopyOldSol(q_temp,q_);  // Copying level 0 solution and saving it

    // Step1:
    //-----------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp1,q_);

    // Step2:
    //------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp1[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp2,q_);

    // Step3:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp2[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp3,q_);

    // Step4:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp3[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp4,q_);

    // Step5:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp4[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp5,q_);

    // Step6:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp5[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp6,q_);

    // Step7:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp6[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);
    CopyOldSol(q_temp7,q_);

    // Step8:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp7[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    // Step9:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
      for(k=0; k<Ndof; k++)
        q_[j][k] = coeff[0]*q_temp[j][k]
                  +coeff[1]*q_temp1[j][k]
                  +coeff[2]*q_temp2[j][k]
                  +coeff[3]*q_temp3[j][k]
                  +coeff[4]*q_temp4[j][k]
                  +coeff[5]*q_temp5[j][k]
                  +coeff[6]*q_temp6[j][k]
                  +coeff[7]*q_temp7[j][k]
                  +coeff[8]*(q_[j][k]+dt_*resid[j][k]);

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::classicRK4(double **q_){

    register int j;

    int k;

    CopyOldSol(q_temp,q_);  // Copying level 0 solution and saving it
    CopyOldResid(resid_temp0,resid);  // Copying level 0 residual and saving it

    // Step1:
    //-----------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] + 0.5 * dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    CopyOldResid(resid_temp1,resid);  // Copying level 1 residual and saving it

    // Step2:
    //------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] + 0.5 * dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    CopyOldResid(resid_temp2,resid);  // Copying level 2 residual and saving it

    // Step3:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    // Step4:
    //--------------
    for(j=0; j<simdata->Nelem_; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k]
                    + (dt_/6.) * ( resid_temp0[j][k] + 2*resid_temp1[j][k]
                                   + 2*resid_temp2[j][k] + resid[j][k] );

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::ComputeInitialResid(double **qn_){

    space_solver->UpdateResid(resid,qn_);

    return;
}

void ExplicitTimeSolver::Reset_time_solver(){

    emptyarray(simdata->Nelem_,resid);
    emptyarray(simdata->Nelem_,q_temp);
    emptyarray(simdata->Nelem_,q_temp1);
    emptyarray(simdata->Nelem_,q_temp2);
    emptyarray(simdata->Nelem_,q_temp3);
    emptyarray(simdata->Nelem_,q_temp4);
    emptyarray(simdata->Nelem_,q_temp5);
    emptyarray(simdata->Nelem_,q_temp6);
    emptyarray(simdata->Nelem_,q_temp7);

    emptyarray(simdata->Nelem_,resid_temp0);
    emptyarray(simdata->Nelem_,resid_temp1);
    emptyarray(simdata->Nelem_,resid_temp2);

    return;
}

