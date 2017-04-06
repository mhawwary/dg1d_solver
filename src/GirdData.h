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

        return;
    }

    void print_grid(){

        register int i;

        cout << "\nPrint grid: "<<endl;
        for(i=0; i<Nfaces; i++)
            _print(i,X[i]);

        return;
    }


    void Reset_(){

        emptyarray(X);
        emptyarray(h_j);

        return;
    }

};

//void GridData::set_grid_param(const SimData& simdata_){

//    Nelem = simdata_.Nelem_;

//    x0 = simdata_.x0_;
//    xf = simdata_.xf_;

//    uniform = simdata_.uniform_;

//    refine_level=simdata_.refine_level_;

//    Nfaces = Nelem+1;

//    X = new double[Nfaces];

//    if(uniform==1) dx = (xf-x0)/Nelem;

//    h_j = new double [Nelem];

//    return;
//}

//void GridData::generate_grid(){

//    register int i;

//    for(i=0; i<Nfaces; i++){

//        X[i]   = dx * (i)  + x0 ;  // node 0, element i
//    }

//    for (i=0; i<Nelem; i++){

//        h_j[i]= X[i+1]-X[i];

//        Xc[i] = 0.5 * ( X[i+1]+X[i] );
//    }

//    return;
//}

//void GridData::Reset_(){

//    emptyarray(X);
//    emptyarray(h_j);

//    return;
//}

//void GridData::print_grid(){

//    register int i;

//    cout << "\nPrint grid: "<<endl;
//    for(i=0; i<Nfaces; i++)
//        _print(i,X[i]);

//    return;
//}

#endif
