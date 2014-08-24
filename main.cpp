#include "fdtd.h"
using namespace std;

const double pi = 3.14159265358979, eps0=8.8541878176*1.0E-12, mu0=1.2566370614*1E-6;

int NUM_OF_PROCESS = 4, NUM_OF_THREADS = 4;

double dx, dy, dz, dt;
int nx, ny, nz, nxprop, nyprop, nzprop;
int nxSize,nxStep;
int nt, nt_of_src;
int output_step_t_of_wavefield, output_step_x_of_wavefield, output_step_of_slice;

int nxPML, nyPML, nzPML;
int i, j, k, ii, it;
int m, kapxmax, kapymax, kapzmax;
int nxslice, nyslice, nzslice;
double alpha;
int nsrc, nrec;
int Total_nPML;

double C, ep_r, sigxmax, sigymax, sigzmax;
double CA,CB,CP,CQ,DH,DE,Bx,By,Bz,Ax,Ay,Az;
double Exdiffy, Exdiffz, Eydiffx, Eydiffz, Ezdiffx, Ezdiffy, Hxdiffy, Hxdiffz, Hydiffx, Hydiffz, Hzdiffx, Hzdiffy;

int *slicex, *slicey, *slicez;
Point *src, *rec;
PSI_E *psi_E;
PSI_H *psi_H;

double *sig, *eps, *mu;
double *src_pulse, *gather;

double *Ex, *Ey, *Ez, *psi_Exy, *psi_Exz, *psi_Eyx, *psi_Eyz, *psi_Ezx, *psi_Ezy;
double *Hx, *Hy, *Hz, *psi_Hxy, *psi_Hxz, *psi_Hyx, *psi_Hyz, *psi_Hzx, *psi_Hzy;

double sigx, sigy, sigz, kapx, kapy, kapz;

char rec_kind_name[10];
char src_kind_name[10];

int myRank;
int order = 2;
int step_x, pml_start=-1, pml_end=-1, main_start = -1, main_end = -1;
int *rangex;
vector<int> mysrc, myrec;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int namelen;
int sendcount,recvcount;
int *displs,*scounts;
MPI_Status status;
double *sig_total, *eps_total, *mu_total;
double *E_total,*H_total;
void add_source();
void record_gather();
int EHtype;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &NUM_OF_PROCESS);
    MPI_Get_processor_name(processor_name, &namelen);
    printf( "%s: Hello world from process %d of %d\n", processor_name, myRank, NUM_OF_PROCESS);

	double t1, t2;
	t1 = MPI_Wtime();

    getParameters();

    if (nx % NUM_OF_PROCESS == 0){
        nxSize = nx/NUM_OF_PROCESS+2*order;
        step_x = nx/NUM_OF_PROCESS;
    }
    else {
        if(myRank == 0) printf("Wanning: Can't be divided even!\n");
        nxSize = nx/NUM_OF_PROCESS+2*order + 1;
        step_x = nx/NUM_OF_PROCESS + 1;
    }

    init();
    readSource();

    readReceive();
    readSlice();

    if (myRank == 0){
        sig_total = new double[(nx+2*order) * ny * nz];
        eps_total = new double[(nx+2*order) * ny * nz];
        mu_total = new double[(nx+2*order) * ny * nz];
        E_total = new double[(nx+2*order) * ny * nz];
        //H_total = new double[(nx+2*order) * ny * nz];
        readSig();
        readEps();
        readMu();
    }



    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(eps_total,scounts,displs,MPI_DOUBLE,eps,scounts[myRank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatterv(mu_total,scounts,displs,MPI_DOUBLE,mu,scounts[myRank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatterv(sig_total,scounts,displs,MPI_DOUBLE,sig,scounts[myRank],MPI_DOUBLE,0,MPI_COMM_WORLD);

	double it1, it2, wt1, wt2, it_total = 0, wt_total = 0;

	for (it=0; it<nt; ++it)
	{
		it1 = MPI_Wtime();

        updateEH();

		it2 = MPI_Wtime();

		it_total += it2 - it1;

        if(myRank == 0){
            cout << "Iteration " << it << ", time: " << it2-it1 << endl;
        }

        if (it < nt_of_src){
            add_source();
		}
        record_gather();

		if ( it % output_step_of_slice == 0 && it != 0) 		// Check point
		{

            //MPI_Gatherv(Ey,scounts[myRank],MPI_DOUBLE,E_total,scounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gatherv(rec[0].component,scounts[myRank],MPI_DOUBLE,E_total,scounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);

			wt1 = MPI_Wtime();

			if (myRank == 0){
                output_slice(it/output_step_of_slice,E_total);
			}

			wt2 = MPI_Wtime();
			wt_total += wt2 - wt1;

		}

		if (it % output_step_t_of_wavefield == 0 && it != 0)
        {
            //MPI_Gatherv(Hx,scounts[myRank],MPI_DOUBLE,E_total,scounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gatherv(rec[mysrc[0]].component,scounts[myRank],MPI_DOUBLE,E_total,scounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);

            wt1 = MPI_Wtime();

            if (myRank == 0){
                output_wavefield(it/output_step_t_of_wavefield);
            }

			wt2 = MPI_Wtime();
			wt_total += wt2 - wt1;

        }
	}

	wt1 = MPI_Wtime();

	output_gather();

	wt2 = MPI_Wtime();

	wt_total += wt2 - wt1;

	t2 = MPI_Wtime();

    if(myRank == 0){
        cout << "Write Total Time: " << wt_total << endl;
        cout << "Iteration Total time: " << it_total << endl;
        cout << "Total Time: " << t2 - t1 << endl;
    }

    MPI_Finalize();

    return 0;
}

void add_source(){

    for (ii=0; ii<(int)mysrc.size(); ++ii)
    {
        i = src[ mysrc[ii] ].x;
        j = src[ mysrc[ii] ].y;
        k = src[ mysrc[ii] ].z;
        int i0 = i % step_x + order;
        src[mysrc[ii]].component(i0,j,k) = src[mysrc[ii]].component(i0,j,k) + src_pulse(ii,it);

        //cout << "src_pulse: " << ii <<" "<< it <<" "<< src_pulse(ii,it) << endl;
        //if(ii == 0) cout << "Ey: " << Ey(i0,j,k) << endl;
    }
    MPI_EH(EHtype);

}

void record_gather(){

    for (ii=0; ii<(int)myrec.size(); ++ii)
    {
        i = rec[ myrec[ii] ].x;
        j = rec[ myrec[ii] ].y;
        k = rec[ myrec[ii] ].z;
        int i0 = i % step_x + order;
        gather(ii,it) = rec[myrec[ii]].component(i0,j,k);

    }
}
