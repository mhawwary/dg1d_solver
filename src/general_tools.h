#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <new>
#include <iomanip>
#include<map>
#ifdef _WIN32
 #include <direct.h>       // those are for mkdir and chdir and other directory and filesystem commands
#elif defined __linux__
 #include <sys/types.h>    // those are for mkdir and chdir and other directory and filesystem commands
 #include <sys/stat.h>     // those are for mkdir and chdir and other directory and filesystem commands
#endif

#include<chrono>

#include <unistd.h>

#include <cstdint>

#include<cstring>

#include <stdio.h>

#include"../include/error.h"

//#include<vector>

//#include"omp.h"
//#include<mpi.h>

using namespace std;


template<typename ptr_1D>
void emptyarray(ptr_1D*& A);

template<typename ptr_2D>
void emptyarray(const int rowsize, ptr_2D**& A);

template<typename ptr_3D>
void emptyarray(const int rowsize, const int colsize, ptr_3D***& A);

template<typename T>
void emptypointer(T*& TT);

// Memory freeing functions:

template<typename T>         // free a one element pointer
void emptypointer(T*& TT){

    if(TT!=nullptr) { delete TT; TT=nullptr; }
}

template<typename ptr_1D>    // free 1D pointer
void emptyarray(ptr_1D*& A)
{
    if(A!=nullptr) { delete [] A; A=nullptr; }
}

template<typename ptr_2D>                             // free 2D pointer
void emptyarray(const int rowsize, ptr_2D**& A)
{
    register int i;
    if(A!=nullptr)
    {
        for(i=0;i<rowsize;i++)  emptyarray(A[i]);
        emptyarray(A);
    }

    return;
}

template<typename ptr_3D>                               // free 3D pointer
void emptyarray(const int rowsize, const int colsize, ptr_3D***& A)
{

    register int i;

    if(A!=nullptr)
    {
        for(i=0;i<rowsize;i++)  emptyarray(colsize,A[i]);
        emptyarray(A);
    }

    return;
}


