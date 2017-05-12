#ifndef GRIDDATA_H
#define GRIDDATA_H

#include "general_tools.h"
#include"SimData.hpp"

//struct SimData;

struct GridData {

    double x0=0;
    double xf=1.0;
    double dx=1.0;

    double *X=nullptr;    // mesh X coord
    double *Xc=nullptr;   // mesh X centroid
    double *h_j=nullptr;  // mesh width for each elements

    double *x_exact_ppts=nullptr;
    double *xi_disc=nullptr;
    int N_exact_ppts=100;
    int N_xi_disc_ppts=1;

    int Nelem=1;
    int Nfaces=2;

    int uniform=1;  // 0: for nonuniform mesh elements

    int refine_level=0; // 0: no refinement


    void set_grid_param(const SimData& simdata_){

        Nelem = simdata_.Nelem_;

        x0 = simdata_.x0_;
        xf = simdata_.xf_;

        uniform = simdata_.uniform_;

        refine_level=simdata_.refine_level_;

        Nfaces = Nelem+1;

        X = new double[Nfaces];

        if(uniform==1) dx = (xf-x0)/Nelem;

        h_j = new double [Nelem];

        Xc = new double [Nelem];

        if(simdata_.poly_order_==0){
            N_xi_disc_ppts =1;
        }else if(simdata_.poly_order_==1){
            N_xi_disc_ppts = 2;
        }else if(simdata_.poly_order_>1){
            N_xi_disc_ppts = simdata_.Npplot;
        }

        N_exact_ppts= 150;

        x_exact_ppts = new double[N_exact_ppts];
        xi_disc = new double[N_xi_disc_ppts];

        return;
    }

    void generate_grid(){

        register int i;

        for(i=0; i<Nfaces; i++){

            X[i]   = dx * (i)  + x0 ;  // node 0, element i
        }

        for (i=0; i<Nelem; i++){

            h_j[i]= X[i+1]-X[i];

            Xc[i] = 0.5 * ( X[i+1]+X[i] );
        }

        // Discontinuous per element sampling
        // for plotting a smooth numerical solution:
        double dxi_=2.0/(N_xi_disc_ppts-1);

        for(i=0; i<N_xi_disc_ppts; i++){

            xi_disc[i] = dxi_ * (i) + -1.0 ;
        }

        // New sampling for plotting a smooth exact solution

        double dxx_ = (xf - x0) / (N_exact_ppts-1);

        for(i=0; i<N_exact_ppts; i++){

            x_exact_ppts[i]   = dxx_ * (i)  + x0 ;
        }

        return;
    }

    void Reset_(){

        emptyarray(X);
        emptyarray(h_j);
        emptyarray(Xc);
        emptyarray(x_exact_ppts);
        emptyarray(xi_disc);

        return;
    }

};


#endif
