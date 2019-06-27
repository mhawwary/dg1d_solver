#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <new>
#include <iomanip>
#include<map>
//#include <limits>
#include <cfloat>
#ifdef _WIN32
 #include <direct.h>
#endif
#include <sys/types.h>    // those are for mkdir and chdir and other directory and filesystem commands
#include <sys/stat.h>     // those are for mkdir and chdir and other directory and filesystem commands

#include<chrono>

#include <unistd.h>

#include <cstdint>

#include<cstring>

#include <stdio.h>

#include <time.h>

#include <omp.h>

#include"../include/error.h"

//#include<vector>

#include"omp.h"
//#include<mpi.h>

using namespace std;


void QuickSort(double*& szArray,double*& ssArray , int nLower, int nUpper);
int Partition(double*& szArray, double*& ssArray,int nLower, int nUpper);

void QuickSort3(double*& mainArray,double*& A1 , int*& A2 , int nLower, int nUpper);
int Partition3(double*& mainArray,double*& A1 , int*& A2, int nLower, int nUpper);

void open_inputfile_forreading(std::string& o_readfname_, std::ifstream& o_input_);
std::string remove_extension(const std::string& filename);
std::string remove_from_end_up_to_string(const std::string field,const std::string& filename);
std::string GetFileExtension(const std::string& FileName);
std::string GetFileDirectory(const std::string& FileName);
void line2intdata(const std::string line, std::vector<int> &data, int& n_data_oneline);
void line2doubledata(const std::string line, std::vector<double> &data, int& n_data_oneline);
bool is_a_comment_line(const std::string& line_in);
bool is_a_text_line(const std::string& line_in);
void copyFile(const std::string& fileNameFrom, const std::string& fileNameTo);

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


