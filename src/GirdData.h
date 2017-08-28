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
    double *xi_uniform=nullptr;
    double *X_dump=nullptr; // Equally Spaced mesh nodes for post processing
    int *x_dump_to_elem=nullptr;

    int Nelem=1;
    int Nfaces=2;

    int N_exact_ppts=100;
    int N_xi_disc_ppts=1;
    int N_uniform_pts=2*Nelem;
    int n_uniform_pts_per_elem=2;

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

        n_uniform_pts_per_elem = simdata_.N_uniform_pts_per_elem_;
        N_uniform_pts = (n_uniform_pts_per_elem-1) * Nelem ;
        X_dump = new double[N_uniform_pts];
        xi_uniform = new double[n_uniform_pts_per_elem];

        if(uniform==1) dx = (xf-x0)/Nelem;

        h_j = new double [Nelem];

        Xc = new double [Nelem];

        if(simdata_.poly_order_==0){
            N_xi_disc_ppts =2;
        }else if(simdata_.poly_order_==1){
            N_xi_disc_ppts = 2;
        }else if(simdata_.poly_order_>1){
            N_xi_disc_ppts = simdata_.Npplot;
        }

        N_exact_ppts= simdata_.N_exact_plot_pts+1;

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

        if(N_xi_disc_ppts==1){

            xi_disc[0] = 0.0;

        }else{

            double dxi_=2.0/(N_xi_disc_ppts-1);

            for(i=0; i<N_xi_disc_ppts; i++)
                xi_disc[i] = dxi_ * (i) + -1.0 ;
        }

        // Globally Equally spaced points for plottings:

        double dxi_u=2.0/(n_uniform_pts_per_elem-1);

        for(i=0; i<n_uniform_pts_per_elem; i++)
            xi_uniform[i] = dxi_u * (i) + -1.0 ;

//        double dx_eq_ = (xf-x0)/(N_equal_spaced);
//        int eID_=0;

//        for(i=0; i<N_equal_spaced+1; i++){
//            X_dump[i] = dx_eq_ *(i) + x0;

//            if(X_dump[i]>X[eID_+1]){
//                eID_++;
//                x_dump_to_elem[i] = eID_;

//            }else if(X_dump[i]==X[eID_+1]){
//                eID_++;
//                x_dump_to_elem[i] = -100;

//            }else{
//                x_dump_to_elem[i] = eID_;
//            }
//        }

//        x_dump_to_elem[N_equal_spaced]=Nelem-1;

        // New sampling for plotting a smooth exact solution

        double dxx_ = (xf - x0) / (N_exact_ppts-1);

        for(i=0; i<N_exact_ppts; i++){

            x_exact_ppts[i]   = dxx_ * (i)  + x0 ;
        }

        return;
    }

    void Reset_(){

        emptyarray(X);
        emptyarray(X_dump);
        emptyarray(x_dump_to_elem);
        emptyarray(h_j);
        emptyarray(Xc);
        emptyarray(x_exact_ppts);
        emptyarray(xi_disc);
        emptyarray(xi_uniform);

        return;
    }

};


#endif
