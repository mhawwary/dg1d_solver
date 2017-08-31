#ifndef QUADRATURE_H
#define QUADRATURE_H

#include"general_tools.h"

struct GaussQuad {

public:

    int Nq=1;  // order of quadrature which is accurate for 2(Nq-1)

    double *Gaus_pts=nullptr;  // Gauss eta=xi Points
    //double *Gaus_X=nullptr;    // Gauss X coord
    double *Gaus_wts=nullptr;  // Gauss weights

public:
    //virtual Face& operator =(const Face& RFace)
    //GausQuad(const int& Nq_);
    //~GausQuad(void);

    void setup_quadrature(const int &Nq_){

        Nq = Nq_;

        Gaus_pts = new double[Nq+1];
        Gaus_wts = new double[Nq+1];

        //Gaus_X  = new double[Nq+1];

        define_gauss_quadrature();

        return;

    }

    void Reset_quad(){

        emptyarray(Gaus_pts);
        emptyarray(Gaus_wts);

        return;
    }

protected:
    void define_gauss_quadrature(){

        if(Nq==1){
            Gaus_pts[0] = 0.0 ;

            Gaus_wts[0] = 2.0 ;

        }else if(Nq==2){
            Gaus_pts[0] = -sqrt(1.0/3) ;
            Gaus_pts[1] =  sqrt(1.0/3) ;

            Gaus_wts[0] =  1.0 ;
            Gaus_wts[1] =  1.0 ;

        }else if(Nq==3){
            Gaus_pts[0] = -sqrt(3.0/5) ;
            Gaus_pts[1] = 0.0  ;
            Gaus_pts[2] = sqrt(3.0/5) ;

            Gaus_wts[0] =  5.0/9 ;
            Gaus_wts[1] =  8.0/9 ;
            Gaus_wts[2] =  5.0/9 ;

        }else if(Nq==4){
            Gaus_pts[0] =  -sqrt(3.0/7+2.0/7*sqrt(6.0/5)) ;
            Gaus_pts[1] =  -sqrt(3.0/7-2.0/7*sqrt(6.0/5)) ;
            Gaus_pts[2] =   sqrt(3.0/7-2.0/7*sqrt(6.0/5)) ;
            Gaus_pts[3] =   sqrt(3.0/7+2.0/7*sqrt(6.0/5)) ;

            Gaus_wts[0] =  (18-sqrt(30))/36.0  ;
            Gaus_wts[1] =  (18+sqrt(30))/36.0  ;
            Gaus_wts[2] =  (18+sqrt(30))/36.0  ;
            Gaus_wts[3] =  (18-sqrt(30))/36.0  ;

        }else if(Nq==5){
            Gaus_pts[0] =  -(1.0/3)*sqrt(5+2*sqrt(10.0/7)) ;
            Gaus_pts[1] =  -(1.0/3)*sqrt(5-2*sqrt(10.0/7)) ;
            Gaus_pts[2] =             0              ;
            Gaus_pts[3] =   (1.0/3)*sqrt(5-2*sqrt(10.0/7)) ;
            Gaus_pts[4] =   (1.0/3)*sqrt(5+2*sqrt(10.0/7)) ;

            Gaus_wts[0] =  (322-13*sqrt(70))/900.0 ;
            Gaus_wts[1] =  (322+13*sqrt(70))/900.0 ;
            Gaus_wts[2] =  128.0/225 ;
            Gaus_wts[3] =  (322+13*sqrt(70))/900.0 ;
            Gaus_wts[4] =  (322-13*sqrt(70))/900.0;

        }else if(Nq==6){
            Gaus_pts[0] =  -9.3246951420315202781230155449399e-01L;
            Gaus_pts[1] =  -6.6120938646626451366139959501991e-01L;
            Gaus_pts[2] =  -2.3861918608319690863050172168071e-01L;
            Gaus_pts[3] =  -Gaus_pts[2];
            Gaus_pts[4] =  -Gaus_pts[1];
            Gaus_pts[5] =  -Gaus_pts[0];

            Gaus_wts[0] =  1.7132449237917034504029614217273e-01L;
            Gaus_wts[1] =  3.6076157304813860756983351383772e-01L;
            Gaus_wts[2] =  4.6791393457269104738987034398955e-01L;
            Gaus_wts[3] =  Gaus_wts[2];
            Gaus_wts[4] =  Gaus_wts[1];
            Gaus_wts[5] =  Gaus_wts[0];

        }else if(Nq==7){
            Gaus_pts[0] =  -9.4910791234275852452618968404785e-01L;
            Gaus_pts[1] =  -7.4153118559939443986386477328079e-01L;
            Gaus_pts[2] =  -4.0584515137739716690660641207696e-01L;
            Gaus_pts[3] =   0.0;
            Gaus_pts[4] =  -Gaus_pts[2];
            Gaus_pts[5] =  -Gaus_pts[1];
            Gaus_pts[6] =  -Gaus_pts[0];

            Gaus_wts[0] =  1.2948496616886969327061143267908e-01L;
            Gaus_wts[1] =  2.7970539148927666790146777142378e-01L;
            Gaus_wts[2] =  3.8183005050511894495036977548898e-01L;
            Gaus_wts[3] =  4.1795918367346938775510204081633e-01L;
            Gaus_wts[4] =  Gaus_wts[2];
            Gaus_wts[5] =  Gaus_wts[1];
            Gaus_wts[6] =  Gaus_wts[0];

        }else if(Nq==8){
            Gaus_pts[0] =  -9.6028985649753623168356086856947e-01L;
            Gaus_pts[1] =  -7.9666647741362673959155393647583e-01L;
            Gaus_pts[2] =  -5.2553240991632898581773904918925e-01L;
            Gaus_pts[3] =  -1.8343464249564980493947614236018e-01L;
            Gaus_pts[4] =  -Gaus_pts[3];
            Gaus_pts[5] =  -Gaus_pts[2];
            Gaus_pts[6] =  -Gaus_pts[1];
            Gaus_pts[7] =  -Gaus_pts[0];

            Gaus_wts[0] =  1.0122853629037625915253135430996e-01L;
            Gaus_wts[1] =  2.2238103445337447054435599442624e-01L;
            Gaus_wts[2] =  3.1370664587788728733796220198660e-01L;
            Gaus_wts[3] =  3.6268378337836198296515044927720e-01L;
            Gaus_wts[4] =  Gaus_wts[3];
            Gaus_wts[5] =  Gaus_wts[2];
            Gaus_wts[6] =  Gaus_wts[1];
            Gaus_wts[7] =  Gaus_wts[0];

        }else{

            cout<< "\n Gauss Quadrature with order of:  "
                <<Nq<<"  is not implemented in the quadrature rule file\n";
        }

        return;
    }
};

#endif

