#include "fdtd.h"
using namespace std;
static char src_file[MAX_NAME_LEN], rec_file[MAX_NAME_LEN];
static char eps_file[MAX_NAME_LEN], mu_file[MAX_NAME_LEN], sigma_file[MAX_NAME_LEN], slice_file[MAX_NAME_LEN];

char dummy_str[10];
double dummy;


void getParameters()
{
    char str[MAX_NAME_LEN];
    FILE *fp;
    fp = fopen("./par.in","r");
    //ifstream fin("./par.in");

    if (fp==NULL)
    {
        printf("Open par.in error!");
        exit(1);
    }

    //fin.getline(str, sizeof(str));
    fgets(str,sizeof(str),fp);
    //cout << "str = " << str << endl;

    fgets(str,sizeof(str),fp);
    sscanf(str, "%lf %lf %lf %lf", &dx, &dy, &dz, &dt);

    fgets(str,sizeof(str),fp);
    fgets(str,sizeof(str),fp);
    sscanf(str, "%d %d %d %d", &nx, &ny, &nz, &nt);

    fgets(str,sizeof(str),fp);
    fgets(str,sizeof(str),fp);
    sscanf(str, "%d", &nt_of_src);

    fgets(str,sizeof(str),fp);
    fgets(str,sizeof(str),fp);
    sscanf(str, "%d %d", &output_step_t_of_wavefield, &output_step_x_of_wavefield);

    fgets(str,sizeof(str),fp);
    fgets(str,sizeof(str),fp);
    sscanf(str, "%d", &output_step_of_slice);

    fgets(str,sizeof(str),fp);
    fgets(str,sizeof(str),fp);
    sscanf(str, "%d %d %d", &nxPML, &nyPML, &nzPML);

    fgets(str,sizeof(str),fp);
    fgets(str,sizeof(str),fp);
    sscanf(str, "%d %d %d %d %lf", &m, &kapxmax, &kapymax, &kapzmax, &alpha);

    fgets(str,sizeof(str),fp);
    fscanf(fp,"%s\n",src_file);

    fgets(str,sizeof(str),fp);
    fscanf(fp,"%s\n",rec_file);

    fgets(str,sizeof(str),fp);
    fscanf(fp,"%s\n",eps_file);

    fgets(str,sizeof(str),fp);
    fscanf(fp,"%s\n",mu_file);

    fgets(str,sizeof(str),fp);
    fscanf(fp,"%s\n",sigma_file);

    fgets(str,sizeof(str),fp);
    fscanf(fp,"%s\n",slice_file);

    fclose(fp);
    if(myRank==0){
    printf("%s\n",src_file);
    printf("%s\n",rec_file);
    printf("%s\n",eps_file);
    printf("%s\n",mu_file);
    printf("%s\n",sigma_file);
    }
}

void readSource()
{
    if(myRank == 0) cout << "call read source" << endl;
    char read[MAX_NAME_LEN] = "./";
    FILE *fp;
    fp = fopen(strcat(read, src_file),"r");
    if (fp==NULL)
    {
        printf("Open src.in error!");
        exit(1);
    }

    fscanf(fp,"%d",&nsrc);
    src = new Point[nsrc];

    fscanf(fp,"%d",&nt_of_src);
    src_pulse = new double[nsrc*nt_of_src];

    for (i=0; i<nsrc; ++i){
        fscanf(fp,"%d %d %d %s",&(src[i].x),&(src[i].y),&(src[i].z),dummy_str);
        if (strcmp(dummy_str,"Ex")==0) {src[i].component = Ex; EHtype = 1;}
        if (strcmp(dummy_str,"Ey")==0) {src[i].component = Ey; EHtype = 1;}
        if (strcmp(dummy_str,"Ez")==0) {src[i].component = Ez; EHtype = 1;}
        if (strcmp(dummy_str,"Hx")==0) {src[i].component = Hx; EHtype = 2;}
        if (strcmp(dummy_str,"Hy")==0) {src[i].component = Hy; EHtype = 2;}
        if (strcmp(dummy_str,"Hz")==0) {src[i].component = Hz; EHtype = 2;}
    }

    int ii = 0;

    for (i=0; i<nsrc; ++i){      // get the sources to be processed by this process
        if (src[i].x>=rangex[myRank] && src[i].x<rangex[myRank+1]){

            mysrc.push_back(i);

            for ( j=0; j<nt_of_src; ++j){
                //cout << "***************" << endl;
                fscanf(fp,"%lf",&(src_pulse(ii,j)));
                //cout << scientific << src_pulse(ii,j) << endl;

            }
            ii += 1;
        }
        else{
            for ( j=0; j<nt_of_src; ++j){
                fscanf(fp,"%lf",&dummy);
            }
        }
    }

    //cout << "Number of sources: " << mysrc.size() << " in myRank: " << myRank << endl;
    fclose(fp);

    if(myRank == 0) cout << "call read source end" << endl;
}

void readReceive()
{
    if(myRank == 0) cout << "call read recv" << endl;
    char read[MAX_NAME_LEN] = "./";
    FILE *fp;
    fp = fopen(strcat(read, rec_file),"r");
    if (fp==NULL)
    {
        printf("Open rec.in error!");
        exit(1);
    }

    fscanf(fp,"%d",&nrec);
    rec = new Point[nrec];

    gather = new double[nrec*nt];

    for (i=0; i<nrec; ++i){
        fscanf(fp,"%d %d %d %s",&(rec[i].x),&(rec[i].y),&(rec[i].z),dummy_str);
        if (strcmp(dummy_str,"Ex")==0) rec[i].component = Ex;
        if (strcmp(dummy_str,"Ey")==0) rec[i].component = Ey;
        if (strcmp(dummy_str,"Ez")==0) rec[i].component = Ez;
        if (strcmp(dummy_str,"Hx")==0) rec[i].component = Hx;
        if (strcmp(dummy_str,"Hy")==0) rec[i].component = Hy;
        if (strcmp(dummy_str,"Hz")==0) rec[i].component = Hz;
    }

    for (i=0; i<nrec; ++i)      // get the receive points to be processed by this process
        if (rec[i].x>=rangex[myRank] && rec[i].x<rangex[myRank+1])
            myrec.push_back(i);
    //cout << "Number of receivers: "<< myrec.size() << " in myRank: " << myRank << endl;

    fclose(fp);
    if(myRank == 0) cout << "call read recv end" << endl;
}

void readSig()
{
    if(myRank == 0) cout << "call read sig" << endl;
    char read[MAX_NAME_LEN] = "./";
    FILE *fp;
    fp = fopen(strcat(read, sigma_file),"r");
    if (fp==NULL)
    {
        printf("Open sig.in error!");
        exit(1);
    }

    for( i=0; i<nx; ++i)
        for( j=0; j<ny; ++j)
            for( k=0; k<nz; ++k){
                fscanf(fp,"%lf",&(sig_total(i+order,j,k)));
                }
    fclose(fp);
}

void readEps()
{
    if(myRank == 0) cout << "call read eps" << endl;
    char read[MAX_NAME_LEN] = "./";
    FILE *fp;
    fp = fopen(strcat(read, eps_file),"r");
    if (fp==NULL)
    {
        printf("Open eps.in error!");
        exit(1);
    }
    for( i=0; i<nx; ++i)
        for( j=0; j<ny; ++j)
            for( k=0; k<nz; ++k){
		fscanf(fp,"%lf",&(dummy));
                eps_total(i+order,j,k) = dummy * eps0;
                }
    fclose(fp);
}

void readMu()
{
    if(myRank == 0) cout << "call read mu" << endl;
    char read[MAX_NAME_LEN] = "./";
    FILE *fp;
    fp = fopen(strcat(read, mu_file),"r");
    if (fp==NULL)
    {
        printf("Open mu.in error!");
        exit(1);
    }
    for( i=0; i<nx; ++i)
        for( j=0; j<ny; ++j)
            for( k=0; k<nz; ++k){
                fscanf(fp,"%lf",&(dummy));
                mu_total(i+order,j,k) = dummy * mu0;
                }
    fclose(fp);
}

void readSlice()
{
    if(myRank == 0) cout << "call read slice" << endl;
    char read[MAX_NAME_LEN] = "./";
    FILE *fp;
    fp = fopen(strcat(read, slice_file),"r");
    if (fp==NULL)
    {
        printf("Open slice.in error!");
        exit(1);
    }

    fscanf(fp,"%d %d %d",&nxslice,&nyslice,&nzslice);
    slicex = new int[nxslice];
    slicey = new int[nyslice];
    slicez = new int[nzslice];


    for ( i=0; i<nxslice; ++i)
	fscanf(fp,"%d %s\n",&(slicex[i]),dummy_str);
    for ( i=0; i<nyslice; ++i)
	fscanf(fp,"%d %s\n",&(slicey[i]),dummy_str);
    for ( i=0; i<nzslice; ++i)
	fscanf(fp,"%d %s\n",&(slicez[i]),dummy_str);

    fclose(fp);
}


void init_psi()
{
    memset(psi_E,0,sizeof(psi_E));
    memset(psi_H,0,sizeof(psi_H));
    int dum = 0;
    for (i=0;i<nx;i++)
    {
        for (j=0;j<ny;j++)
        {
            for (k=0;k<nz;k++)
            {
                if (i>=nxPML && i<nx-nxPML && j>=nyPML && j<ny-nyPML && k>=nzPML && k<nz-nzPML) continue;
                if ( pml_start==-1 && i>=rangex[myRank] )
                    pml_start = dum;
                if ( pml_start!=-1 && pml_end==-1 && i>=rangex[myRank+1] )
                    pml_end = dum;
                psi_E[dum].x = i;
                psi_E[dum].y = j;
                psi_E[dum].z = k;

                psi_H[dum].x = i;
                psi_H[dum].y = j;
                psi_H[dum].z = k;

                dum = dum + 1;


            }
        }
    }
    if (pml_end==-1)        // the last process
        pml_end = dum;
    //if (myRank == NUM_OF_PROCESS-1 ) pml_start = pml_start -10;

    cout << myRank << " PML: " << pml_start << " " << pml_end << endl; 
}

void init()
{
    if(myRank == 0) cout << "call init" << endl;

    sendcount= order*ny*nz;
    recvcount = sendcount;

    rangex = new int[NUM_OF_PROCESS+1];
    scounts = new int[NUM_OF_PROCESS];
    displs = new int[NUM_OF_PROCESS];

    rangex[0] = 0;
    displs[0] = 0;

    for(i=1;i<NUM_OF_PROCESS;i++){
        if(i<=nx%NUM_OF_PROCESS){
            rangex[i] = rangex[i-1] + (nx/NUM_OF_PROCESS+1);
            scounts[i-1] = ((nx/NUM_OF_PROCESS+1) + 2*order)*ny*nz;
            displs[i] = (rangex[i])*ny*nz;
        }
        else{
            rangex[i] = rangex[i-1] + nx/NUM_OF_PROCESS;
            scounts[i-1] = ((nx/NUM_OF_PROCESS) + 2*order)*ny*nz;
            displs[i] = (rangex[i])*ny*nz;
        }
    }
    rangex[NUM_OF_PROCESS] = nx;
    scounts[NUM_OF_PROCESS-1] = ((nx/NUM_OF_PROCESS) + 2*order)*ny*nz;

    step_x = rangex[myRank+1] - rangex[myRank];
    nxSize = step_x+2*order;
    //if(myRank==NUM_OF_PROCESS-1) step_x = step_x +1;


    if (nxPML<=rangex[myRank])	main_start = order;
    else if (nxPML>=rangex[myRank+1])	main_start=INF;
    else main_start = nxPML - rangex[myRank] + order;

    if (rangex[myRank+1] <= nx-nxPML) main_end = order+step_x;
    else if (rangex[myRank] >= nx-nxPML) main_end=-INF;
    else main_end = order + step_x - (rangex[myRank+1] - (nx-nxPML));


    printf("******:%d step_x = %d, rangex[%d] = %d, main_start=%d, main_end=%d scounts[%d]=%d\n", \
        myRank,step_x, myRank, rangex[myRank],main_start,main_end,myRank,scounts[myRank]/ny/nz);


    nxprop = 2*nx-1;
    nyprop = 2*ny-1;
    nzprop = 2*nz-1;

    int EHFieldSize = nxSize * ny * nz * sizeof(double);
    Total_nPML = nx*ny*nz - (nx-2*nxPML)*(ny-2*nyPML)*(nz-2*nzPML);

    psi_E = new PSI_E[Total_nPML];
    psi_H = new PSI_H[Total_nPML];


    sig = new double[nxSize * ny * nz];
    eps = new double[nxSize * ny * nz];
    mu = new double[nxSize * ny * nz];
    //cout <<"nxSize: " <<nxSize << endl;

    for ( i=0; i<nxSize; i++)
        for ( j=0; j<ny; j++)
            for ( k=0; k<nz; k++)
            {
                eps(i,j,k) = 9.0 * eps0;
                mu(i,j,k) = mu0;
                sig(i,j,k) = 0.0;
            }

    Ex = new double[nxSize*ny*nz];
    Ey = new double[nxSize*ny*nz];
    Ez = new double[nxSize*ny*nz];
    Hx = new double[nxSize*ny*nz];
    Hy = new double[nxSize*ny*nz];
    Hz = new double[nxSize*ny*nz];
    memset(Ex, 0, EHFieldSize);
    memset(Ey, 0, EHFieldSize);
    memset(Ez, 0, EHFieldSize);
    memset(Hx, 0, EHFieldSize);
    memset(Hy, 0, EHFieldSize);
    memset(Hz, 0, EHFieldSize);

    Exdiffy = 0.0;
    Exdiffz = 0.0;
    Eydiffx = 0.0;
    Eydiffz = 0.0;
    Ezdiffx = 0.0;
    Ezdiffy = 0.0;
    Hxdiffy = 0.0;
    Hxdiffz = 0.0;
    Hydiffx = 0.0;
    Hydiffz = 0.0;
    Hzdiffx = 0.0;
    Hzdiffy = 0.0;
    CA = 0.0;
    CB = 0.0;
    CP = 0.0;
    CQ = 0.0;
    DH = 0.0;
    DE = 0.0;
    Bx = 0.0;
    By = 0.0;
    Bz = 0.0;
    Ax = 0.0;
    Ay = 0.0;
    Az = 0.0;
    init_psi();
    if(myRank == 0) cout << "call init end" << endl;
}


