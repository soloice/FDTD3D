#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "omp.h"
#include <ctime>
#include <mpi.h>
#include <vector>


#define MAX_NUM_OF_PROCESS 64
#define MAX_NAME_LEN 300

#define Hx(i,j,k) Hx[( (i)*ny + j )*nz + k]
#define Hy(i,j,k) Hy[( (i)*ny + j )*nz + k]
#define Hz(i,j,k) Hz[( (i)*ny + j )*nz + k]

#define Ex(i,j,k) Ex[( (i)*ny + j )*nz + k]
#define Ey(i,j,k) Ey[( (i)*ny + j )*nz + k]
#define Ez(i,j,k) Ez[( (i)*ny + j )*nz + k]

#define E_total(i,j,k) E_total[( (i)*ny + j )*nz + k]
#define H_total(i,j,k) H_total[( (i)*ny + j )*nz + k]
#define EH_total(i,j,k) EH_total[( (i)*ny + j )*nz + k]


#define sig(i,j,k) sig[( (i)*ny + j )*nz + k]
#define eps(i,j,k) eps[( (i)*ny + j )*nz + k]
#define mu(i,j,k) mu[( (i)*ny + j )*nz + k]
#define sig_total(i,j,k) sig_total[( (i)*ny + j )*nz + k]
#define eps_total(i,j,k) eps_total[( (i)*ny + j )*nz + k]
#define mu_total(i,j,k) mu_total[( (i)*ny + j )*nz + k]

#define src_pulse(i,j) src_pulse[(i)*nt_of_src + j]
#define gather(i,j) gather[(i)*nt + j]
#define component(i,j,k) component[( (i)*ny + j )*nz + k]
#define INF  100000

using std::vector;

extern const double pi, eps0, mu0;

extern double dx, dy, dz, dt;
extern int nx, ny, nz, nxprop, nyprop, nzprop;
extern int nxSize,step_x;
extern int nt, nt_of_src;
extern int output_step_t_of_wavefield, output_step_x_of_wavefield, output_step_of_slice;
extern int nxPML, nyPML, nzPML;
extern int i, j, k, ii, it;
extern int m, kapxmax, kapymax, kapzmax;
extern double alpha;
extern int nsrc, nrec;
extern int nxslice, nyslice, nzslice;

extern double sigxmax, sigymax, sigzmax;
extern double sigx, sigy, sigz, kapx, kapy, kapz;
extern double CA,CB,CP,CQ,DH,DE,Bx,By,Bz,Ax,Ay,Az;
extern double Exdiffy, Exdiffz, Eydiffx, Eydiffz, Ezdiffx, Ezdiffy, Hxdiffy, Hxdiffz, Hydiffx, Hydiffz, Hzdiffx, Hzdiffy;

typedef struct Point{
    int x, y, z;
    //char component[10];
    double *component;
} Point;

extern int Total_nPML;

typedef struct PSI_E{
    int x;
    int y;
    int z;
    double psi_Exy;
    double psi_Eyx;
    double psi_Eyz;
    double psi_Ezy;
    double psi_Ezx;
    double psi_Exz;
} PSI_E;
typedef struct PSI_H{
    int x;
    int y;
    int z;
    double psi_Hxy;
    double psi_Hyx;
    double psi_Hyz;
    double psi_Hzy;
    double psi_Hzx;
    double psi_Hxz;
} PSI_H;


extern PSI_E *psi_E;
extern PSI_H *psi_H;


extern Point *src, *rec;
extern int *slicex,*slicey,*slicez;

extern double *sig, *eps, *mu;
extern double *sig_total, *eps_total, *mu_total;
extern double *src_pulse, *gather;
extern double *Ex, *Ey, *Ez;
extern double *Hx, *Hy, *Hz;
extern double *E_total,*H_total;

extern char rec_kind_name[10];
extern char src_kind_name[10];
extern int NUM_OF_THREADS;

extern int pml_start,pml_end;
extern int *rangex;
extern int myRank,NUM_OF_PROCESS;
extern int step_x, pml_start, pml_end, main_start, main_end;
extern vector<int> mysrc, myrec;
extern int *displs,*scounts;
extern int order;
extern int sendcount,recvcount;
extern MPI_Status status;
extern int EHtype;

void updateEH();
void MPI_EH(int dummy);

void output_gather();
void output_slice(int it,double *EH_total);
void output_wavefield(int it);

void getParameters();
void init();
void readSource();
void readReceive();
void readSig();
void readEps();
void readMu();
void readSlice();
