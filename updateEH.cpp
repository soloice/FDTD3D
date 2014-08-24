#include "fdtd.h"
using namespace std;
double mu_halfx,mu_halfy,mu_halfz;
double sig_halfxy,sig_halfxz,sig_halfyz;
double eps_halfxy,eps_halfxz,eps_halfyz;
//double sigx, sigy, sigz, kapx, kapy, kapz,
double iratio;
int dummy_int;

void updateH_main()
{
    //#pragma omp parallel for num_threads(NUM_OF_THREADS) private(j,k,CQ,mu_halfx,mu_halfy,mu_halfz,\
     Eydiffz,Ezdiffy,Exdiffz,Ezdiffx,Eydiffx,Exdiffy)

    for (i=main_start; i<main_end; ++i)
    {
        for (j=nyPML;j<ny-nyPML;++j)
        {
            for (k=nzPML;k<nz-nzPML;++k)
            {
                mu_halfx = (mu(i,j,k) + mu(i+1,j,k))/2;
                CQ = dt/mu_halfx;

                Eydiffz = 1.0/(24*dz) * ( -Ey(i,j,k+1) + 27*Ey(i,j,k) - 27*Ey(i,j,k-1) + Ey(i,j,k-2));
                Ezdiffy = 1.0/(24*dy) * ( -Ez(i,j+1,k) + 27*Ez(i,j,k) - 27*Ez(i,j-1,k) + Ez(i,j-2,k));
                Hx(i,j,k) = Hx(i,j,k) - CQ*Ezdiffy + CQ*Eydiffz;

                mu_halfy = (mu(i,j,k) + mu(i,j+1,k))/2;
                CQ = dt/mu_halfy;
                Exdiffz = 1.0/(24*dz) * ( -Ex(i,j,k+1) + 27 * Ex(i,j,k) - 27 * Ex(i,j,k-1) + Ex(i,j,k-2));
                Ezdiffx = 1.0/(24*dx) * ( -Ez(i+1,j,k) + 27 * Ez(i,j,k) - 27 * Ez(i-1,j,k) + Ez(i-2,j,k));
                Hy(i,j,k) = Hy(i,j,k) - CQ*Exdiffz + CQ*Ezdiffx;

                mu_halfz = (mu(i,j,k) + mu(i,j,k+1))/2;
                //cout << mu_halfz << endl;
                CQ = dt/mu_halfz;
                //cout << CQ << endl;
                Eydiffx = 1.0/(24*dx)*( -Ey(i+1,j,k) +27*Ey(i,j,k) -27*Ey(i-1,j,k) +Ey(i-2,j,k));
                Exdiffy = 1.0/(24*dy)*( -Ex(i,j+1,k) +27*Ex(i,j,k) -27*Ex(i,j-1,k) +Ex(i,j-2,k));
                Hz(i,j,k) = Hz(i,j,k) - CQ*Eydiffx + CQ*Exdiffy;

            }
        }
    }
}

void updateE_main()
{
    // #pragma omp parallel for num_threads(NUM_OF_THREADS) private(j,k,CA,CB,sig_halfyz,eps_halfyz,\
    sig_halfxz,eps_halfxz,sig_halfxy,eps_halfxy,Hzdiffy,Hydiffz,Hxdiffz,Hzdiffx,Hydiffx,Hxdiffy)

    for (i=main_start; i<main_end; ++i)
    {
        for (j=nyPML;j<ny-nyPML;++j)
        {
            for (k=nzPML;k<nz-nzPML;++k)
            {
                sig_halfyz = (sig(i,j,k) + sig(i,j+1,k+1))/2;
                eps_halfyz = (eps(i,j+1,k+1) + eps(i,j,k))/2;

                CA = (1.0 - sig_halfyz*dt/(2.0*eps_halfyz))/(1.0 + sig_halfyz * dt/(2.0*eps_halfyz));
                CB = (dt/eps_halfyz)/(1.0 + sig_halfyz*dt / (2.0*eps_halfyz));

                Hzdiffy = 1.0/(24*dy)*( -Hz(i,j+2,k) +27*Hz(i,j+1,k) -27*Hz(i,j,k) +Hz(i,j-1,k));
                Hydiffz = 1.0/(24*dz)*( -Hy(i,j,k+2) +27*Hy(i,j,k+1) -27*Hy(i,j,k) +Hy(i,j,k-1));
                Ex(i,j,k) = CA*Ex(i,j,k) + CB*Hzdiffy - CB*Hydiffz;

                sig_halfxz = (sig(i,j,k) + sig(i+1,j,k+1))/2;
                eps_halfxz = (eps(i,j,k) + eps(i+1,j,k+1))/2;
                CA = (1.0 - sig_halfxz*dt/(2.0*eps_halfxz))/(1.0 + sig_halfxz * dt/(2.0*eps_halfxz));
                CB = (dt/eps_halfxz)/(1.0 + sig_halfxz*dt / (2.0*eps_halfxz));

                Hxdiffz = 1.0/(24*dz)*( -Hx(i,j,k+2) +27*Hx(i,j,k+1) -27*Hx(i,j,k) +Hx(i,j,k-1));
                Hzdiffx = 1.0/(24*dx)*( -Hz(i+2,j,k) +27*Hz(i+1,j,k) -27*Hz(i,j,k) +Hz(i-1,j,k));
                Ey(i,j,k) = CA*Ey(i,j,k) + CB*Hxdiffz - CB*Hzdiffx;

                sig_halfxy = (sig(i,j,k) + sig(i+1,j+1,k))/2;
                eps_halfxy = (eps(i,j,k) + eps(i+1,j+1,k))/2;
                CA = (1.0 - sig_halfxy*dt/(2.0*eps_halfxy))/(1.0 + sig_halfxy * dt/(2.0*eps_halfxy));
                CB = (dt/eps_halfxy)/(1.0 + sig_halfxy*dt / (2.0*eps_halfxy));
                Hydiffx = 1.0/(24*dx)*( -Hy(i+2,j,k) +27*Hy(i+1,j,k) -27*Hy(i,j,k) +Hy(i-1,j,k));
                Hxdiffy = 1.0/(24*dy)*( -Hx(i,j+2,k) +27*Hx(i,j+1,k) -27*Hx(i,j,k) +Hx(i,j-1,k));
                Ez(i,j,k) = CA*Ez(i,j,k) + CB*Hydiffx - CB*Hxdiffy;
            }
        }
    }
}


void updateH_PML()
{
    //#pragma omp parallel for num_threads(NUM_OF_THREADS) private(i,j,k,CQ,mu_halfx,mu_halfy,mu_halfz,\
     Eydiffz,Ezdiffy,Exdiffz,Ezdiffx,Eydiffx,Exdiffy,Bx,By,Bz,Ax,Ay,Az,sigx,sigy,sigz,kapx,kapy,kapz,iratio,sigxmax,sigymax,sigzmax)

    for (ii=pml_start; ii<pml_end; ii++)
    {
        i = psi_H[ii].x; j = psi_H[ii].y; k = psi_H[ii].z;
        int i0 = i % step_x + order;
        sigxmax = (m+1) / (150 * pi * sqrt(eps(i0,j,k)/eps0) * dx);
        sigymax = (m+1) / (150 * pi * sqrt(eps(i0,j,k)/eps0) * dy);
        sigzmax = (m+1) / (150 * pi * sqrt(eps(i0,j,k)/eps0) * dz);

        if(i<nxPML){
            iratio = pow( double(2*nxPML-1 - 2*i)/double(2*nxPML-1), m);
            sigx = sigxmax * iratio;
            kapx = 1.0 + (kapxmax - 1.0) * iratio;}
        else if(i>=nx-nxPML){
            iratio = pow( double(2*i-(nxprop-2*nxPML))/double(2*nxPML-1), m );
            sigx = sigxmax * iratio;
            kapx = 1.0 + (kapxmax - 1.0) * iratio;}
        else{
            sigx = 0.0;
            kapx = 1.0;}
        if(j<nyPML){
            iratio = pow( double(2*nyPML-1 - 2*j)/double(2*nyPML-1), m );
            sigy = sigymax * iratio;
            kapy = 1.0 + (kapymax - 1.0) * iratio;}
        else if(j>=ny-nyPML){
            iratio = pow( double(2*j-(nyprop-2*nyPML))/double(2*nyPML-1), m );
            sigy = sigymax * iratio;
            kapy = 1.0 + (kapymax - 1.0) * iratio;}
        else{
            sigy = 0.0;
            kapy = 1.0;}
        if(k<nzPML){
            iratio = pow( double(2*nzPML-1 - 2*k)/double(2*nzPML-1), m );
            sigz = sigzmax * iratio;
            kapz = 1.0 + (kapzmax - 1.0) * iratio;}
        else if(k>=nz-nzPML){
            iratio = pow( double(2*k-(nzprop-2*nzPML))/double(2*nzPML-1), m );
            sigz = sigzmax * iratio;
            kapz = 1.0 + (kapzmax - 1.0) * iratio;}
        else{
            sigz = 0.0;
            kapz = 1.0;}

        if ( j >= 2 && k >= 2 && i <= nx-3 && j <= ny-3 && k <= nz-3 ){
        mu_halfx = (mu(i0,j,k) + mu(i0+1,j,k))/2;
        CQ = dt/mu_halfx;
        Eydiffz = 1.0/(24*dz) * ( -Ey(i0,j,k+1) + 27*Ey(i0,j,k) - 27*Ey(i0,j,k-1) + Ey(i0,j,k-2));
        Ezdiffy = 1.0/(24*dy) * ( -Ez(i0,j+1,k) + 27*Ez(i0,j,k) - 27*Ez(i0,j-1,k) + Ez(i0,j-2,k));
        Hx(i0,j,k) = Hx(i0,j,k) - CQ * Ezdiffy/kapy + CQ*Eydiffz/kapz;

        By = exp( -dt/eps0 * (sigy / kapy + alpha) );
        Bz = exp( -dt/eps0 * (sigz / kapz + alpha) );
        Ay = (sigy / (sigy * kapy + alpha * kapy * kapy + 1E-12)) * (By - 1.0);
        Az = (sigz / (sigz * kapz + alpha * kapz * kapz + 1E-12)) * (Bz - 1.0);
        psi_H[ii].psi_Hxz = Bz*psi_H[ii].psi_Hxz + Az*Eydiffz;
        psi_H[ii].psi_Hxy = By*psi_H[ii].psi_Hxy + Ay*Ezdiffy;
        Hx(i0,j,k) = Hx(i0,j,k) - CQ*psi_H[ii].psi_Hxy + CQ*psi_H[ii].psi_Hxz;
        }
        if ( i >=2 && k >=2 && i <= nx-3 && j <= ny-3 && k <= nz-3){
        mu_halfy = (mu(i0,j,k) + mu(i0,j+1,k))/2;
        CQ = dt/mu_halfy;
        Exdiffz = 1.0/(24*dz) * ( -Ex(i0,j,k+1) + 27 * Ex(i0,j,k) - 27 * Ex(i0,j,k-1) + Ex(i0,j,k-2));
        Ezdiffx = 1.0/(24*dx) * ( -Ez(i0+1,j,k) + 27 * Ez(i0,j,k) - 27 * Ez(i0-1,j,k) + Ez(i0-2,j,k));
        Hy(i0,j,k) = Hy(i0,j,k) - CQ*Exdiffz/kapz + CQ*Ezdiffx/kapx;

        Bx = exp(-dt/eps0 * (sigx/kapx+alpha));
        Bz = exp(-dt/eps0 * (sigz/kapz+alpha));
        Ax = (sigx / (sigx * kapx + alpha * kapx * kapx + 1E-12)) * (Bx - 1.0);
        Az = (sigz / (sigz * kapz + alpha * kapz * kapz + 1E-12)) * (Bz - 1.0);
        psi_H[ii].psi_Hyz = Bz*psi_H[ii].psi_Hyz + Az*Exdiffz;
        psi_H[ii].psi_Hyx = Bx*psi_H[ii].psi_Hyx + Ax*Ezdiffx;
        Hy(i0,j,k) = Hy(i0,j,k) - CQ*psi_H[ii].psi_Hyz + CQ*psi_H[ii].psi_Hyx;
        }
        if ( i>= 2 && j>=2 && i <= nx-3 && j <= ny-3 && k <= nz-3){
        mu_halfz = (mu(i0,j,k) + mu(i0,j,k+1))/2;
        CQ = dt/mu_halfz;
        Eydiffx = 1.0/(24*dx)*( -Ey(i0+1,j,k) +27*Ey(i0,j,k) -27*Ey(i0-1,j,k) +Ey(i0-2,j,k));
        Exdiffy = 1.0/(24*dy)*( -Ex(i0,j+1,k) +27*Ex(i0,j,k) -27*Ex(i0,j-1,k) +Ex(i0,j-2,k));
        Hz(i0,j,k) = Hz(i0,j,k) - CQ*Eydiffx/kapx + CQ*Exdiffy/kapy;

        Bx = exp(-dt/eps0 * (sigx/kapx+alpha));
        By = exp(-dt/eps0 * (sigy/kapy+alpha));
        Ax = (sigx/(sigx*kapx + alpha*kapx*kapx + 1E-12)) * (Bx-1.0);
        Ay = (sigy/(sigy*kapy + alpha*kapy*kapy + 1E-12)) * (By-1.0);
        psi_H[ii].psi_Hzx = Bx*psi_H[ii].psi_Hzx + Ax*Eydiffx;
        psi_H[ii].psi_Hzy = By*psi_H[ii].psi_Hzy + Ay*Exdiffy;
        Hz(i0,j,k) = Hz(i0,j,k) - CQ*psi_H[ii].psi_Hzx + CQ*psi_H[ii].psi_Hzy;
        }
    }
}

void updateE_PML()
{
    //#pragma omp parallel for num_threads(NUM_OF_THREADS) private(i,j,k,CA,CB,sig_halfyz,eps_halfyz,sig_halfxz,eps_halfxz,sig_halfxy,eps_halfxy,\
     Hzdiffy,Hydiffz,Hxdiffz,Hzdiffx,Hydiffx,Hxdiffy,DE,Bx,By,Bz,Ax,Ay,Az,sigx,sigy,sigz,kapx,kapy,kapz,iratio,sigxmax,sigymax,sigzmax)

    for (ii=pml_start; ii<pml_end; ii++)
    {
        i = psi_E[ii].x; j = psi_E[ii].y; k = psi_E[ii].z;

        int i0 = i % step_x + order;

        sigxmax = (m+1) / (150 * pi * sqrt(eps(i0,j,k)/eps0) * dx);
        sigymax = (m+1) / (150 * pi * sqrt(eps(i0,j,k)/eps0) * dy);
        sigzmax = (m+1) / (150 * pi * sqrt(eps(i0,j,k)/eps0) * dz);

        if(i<nxPML){
            iratio = pow( double(2*nxPML-1 - (2*i+1))/double(2*nxPML-1), m);

            sigx = sigxmax * iratio;
            kapx = 1.0 + (kapxmax - 1.0) * iratio;}
        else if(i>=nx-nxPML){
            iratio = pow( double(2*i+1-(nxprop-2*nxPML))/double(2*nxPML-1), m );

            sigx = sigxmax * iratio;
            kapx = 1.0 + (kapxmax - 1.0) * iratio;}
        else{
            sigx = 0.0;
            kapx = 1.0;}
        if(j<nyPML){
            iratio = pow( double(2*nyPML-1 - (2*j+1))/double(2*nyPML-1), m );

            sigy = sigymax * iratio;
            kapy = 1.0 + (kapymax - 1.0) * iratio;}
        else if(j>=ny-nyPML){
            iratio = pow( double(2*j+1-(nyprop-2*nyPML))/double(2*nyPML-1), m );
            sigy = sigymax * iratio;
            kapy = 1.0 + (kapymax - 1.0) * iratio;}
        else{
            sigy = 0.0;
            kapy = 1.0;}

        if(k<nzPML){
            iratio = pow( double(2*nzPML-1 - (2*k+1))/double(2*nzPML-1), m );
            sigz = sigzmax * iratio;
            kapz = 1.0 + (kapzmax - 1.0) * iratio;}
        else if(k>=nz-nzPML){
            iratio = pow( double(2*k+1-(nzprop-2*nzPML))/double(2*nzPML-1), m );
            sigz = sigzmax * iratio;
            kapz = 1.0 + (kapzmax - 1.0) * iratio;}
        else{
            sigz = 0.0;
            kapz = 1.0;}


        if ( j>= 1 && k >= 1 && i <= nx-3 && j <= ny-3 && k <= nz-3){

        sig_halfyz = (sig(i0,j,k) + sig(i0,j+1,k+1))/2;
        eps_halfyz = (eps(i0,j+1,k+1) + eps(i0,j,k))/2;
        CA = (1.0 - sig_halfyz * dt/(2.0*eps_halfyz))/(1.0 + sig_halfyz * dt/(2.0*eps_halfyz));
        CB = (dt/eps_halfyz)/(1.0 + sig_halfyz*dt / (2.0*eps_halfyz));
        Hzdiffy = 1.0/(24*dy)*( -Hz(i0,j+2,k) +27*Hz(i0,j+1,k) -27*Hz(i0,j,k) +Hz(i0,j-1,k));
        Hydiffz = 1.0/(24*dz)*( -Hy(i0,j,k+2) +27*Hy(i0,j,k+1) -27*Hy(i0,j,k) +Hy(i0,j,k-1));
        Ex(i0,j,k) = CA*Ex(i0,j,k) + CB*Hzdiffy/kapy - CB*Hydiffz/kapz;

        DE = dt/eps_halfyz/(1+sig_halfyz*dt/(2.0*eps_halfyz));
        By = exp(-dt/eps0*(sigy/kapy+alpha));
        Bz = exp(-dt/eps0*(sigz/kapz+alpha));
        Ay = (sigy/(sigy*kapy + alpha*kapy*kapy + 1E-12))*(By-1.0);
        Az = (sigz/(sigz*kapz + alpha*kapz*kapz + 1E-12))*(Bz-1.0);
        psi_E[ii].psi_Exy = By*psi_E[ii].psi_Exy + Ay*Hzdiffy;
        psi_E[ii].psi_Exz = Bz*psi_E[ii].psi_Exz + Az*Hydiffz;
        Ex(i0,j,k) = Ex(i0,j,k) + DE*psi_E[ii].psi_Exy - DE*psi_E[ii].psi_Exz;
        }
        if ( i>= 1 && k >= 1 && i <= nx-3 && j <= ny-3 && k <= nz-3){
        sig_halfxz = (sig(i0,j,k) + sig(i0+1,j,k+1))/2;
        eps_halfxz = (eps(i0,j,k) + eps(i0+1,j,k+1))/2;
        CA = (1.0 - sig_halfxz*dt/(2.0*eps_halfxz))/(1.0 + sig_halfxz * dt/(2.0*eps_halfxz));
        CB = (dt/eps_halfxz)/(1.0 + sig_halfxz*dt / (2.0*eps_halfxz));
        Hxdiffz = 1.0/(24*dz)*( -Hx(i0,j,k+2) +27*Hx(i0,j,k+1) -27*Hx(i0,j,k) +Hx(i0,j,k-1));
        Hzdiffx = 1.0/(24*dx)*( -Hz(i0+2,j,k) +27*Hz(i0+1,j,k) -27*Hz(i0,j,k) +Hz(i0-1,j,k));
        Ey(i0,j,k) = CA*Ey(i0,j,k) + CB*Hxdiffz/kapz - CB*Hzdiffx/kapx;
        DE = dt/eps_halfxz/(1+sig_halfxz*dt/(2.0*eps_halfxz));
        Bx = exp(-dt/eps0*(sigx/kapx+alpha));
        Bz = exp(-dt/eps0*(sigz/kapz+alpha));
        Ax = (sigx/(sigx*kapx + alpha*kapx*kapx + 1E-12))*(Bx-1.0);
        Az = (sigz/(sigz*kapz + alpha*kapz*kapz + 1E-12))*(Bz-1.0);
        psi_E[ii].psi_Eyx = Bx*psi_E[ii].psi_Eyx + Ax*Hzdiffx;
        psi_E[ii].psi_Eyz = Bz*psi_E[ii].psi_Eyz + Az*Hxdiffz;
        Ey(i0,j,k) = Ey(i0,j,k) + DE*psi_E[ii].psi_Eyz - DE*psi_E[ii].psi_Eyx;
        }

        if ( i>=1 && j >= 1 && i <= nx-3 && j <= ny-3 && k <= nz-3){
        sig_halfxy = (sig(i0,j,k) + sig(i0+1,j+1,k))/2;
        eps_halfxy = (eps(i0,j,k) + eps(i0+1,j+1,k))/2;
        CA = (1.0 - sig_halfxy*dt/(2.0*eps_halfxy))/(1.0 + sig_halfxy * dt/(2.0*eps_halfxy));
        CB = (dt/eps_halfxy)/(1.0 + sig_halfxy*dt / (2.0*eps_halfxy));
        Hydiffx = 1.0/(24*dx)*( -Hy(i0+2,j,k) +27*Hy(i0+1,j,k) -27*Hy(i0,j,k) +Hy(i0-1,j,k));
        Hxdiffy = 1.0/(24*dy)*( -Hx(i0,j+2,k) +27*Hx(i0,j+1,k) -27*Hx(i0,j,k) +Hx(i0,j-1,k));
        Ez(i0,j,k) = CA*Ez(i0,j,k) + CB*Hydiffx/kapx - CB*Hxdiffy/kapy;

        DE = dt/eps_halfxy/(1+sig_halfxy*dt/(2.0*eps_halfxy));
        Bx = exp(-dt/eps0*(sigx/kapx+alpha));
        By = exp(-dt/eps0*(sigy/kapy+alpha));
        Ax = (sigx/(sigx*kapx + alpha*kapx*kapx + 1E-12))*(Bx-1.0);
        Ay = (sigy/(sigy*kapy + alpha*kapy*kapy + 1E-12))*(By-1.0);
        psi_E[ii].psi_Ezx = Bx*psi_E[ii].psi_Ezx + Ax*Hydiffx;
        psi_E[ii].psi_Ezy = By*psi_E[ii].psi_Ezy + Ay*Hxdiffy;
        Ez(i0,j,k) = Ez(i0,j,k) + DE*psi_E[ii].psi_Ezx - DE*psi_E[ii].psi_Ezy;

        }
    }

}


void MPI_E()
{
    //Passing message to right & receiving message from right
    if (myRank<NUM_OF_PROCESS-1)
    {
        MPI_Sendrecv(&Ey(step_x,0,0), sendcount, MPI_DOUBLE, myRank+1, 1, &Ey(step_x+order,0,0), recvcount, MPI_DOUBLE, myRank+1, 3, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&Ez(step_x,0,0), sendcount, MPI_DOUBLE, myRank+1, 2, &Ez(step_x+order,0,0), recvcount, MPI_DOUBLE, myRank+1, 4, MPI_COMM_WORLD, &status);
    }

    //Passing message to left & receiving message from left
    if (myRank>0)
    {
        MPI_Sendrecv(&Ey(order,0,0), sendcount, MPI_DOUBLE, myRank-1, 3, &Ey(0,0,0), recvcount, MPI_DOUBLE, myRank-1, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&Ez(order,0,0), sendcount, MPI_DOUBLE, myRank-1, 4, &Ez(0,0,0), recvcount, MPI_DOUBLE, myRank-1, 2, MPI_COMM_WORLD, &status);
    }
}

void MPI_H()
{
    //Passing message to right & receiving message from right
    if (myRank<NUM_OF_PROCESS-1)
    {
        //cout << "H right comm" << endl;
        MPI_Sendrecv(&Hy(step_x,0,0), sendcount, MPI_DOUBLE, myRank+1, 1, &Hy(step_x+order,0,0), recvcount, MPI_DOUBLE, myRank+1, 3, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&Hz(step_x,0,0), sendcount, MPI_DOUBLE, myRank+1, 2, &Hz(step_x+order,0,0), recvcount, MPI_DOUBLE, myRank+1, 4, MPI_COMM_WORLD, &status);
    }

    //Passing message to left & receiving message from left
    if (myRank>0)
    {
        //cout << "H left comm" << endl;
        MPI_Sendrecv(&Hy(order,0,0), sendcount, MPI_DOUBLE, myRank-1, 3, &Hy(0,0,0), recvcount, MPI_DOUBLE, myRank-1, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&Hz(order,0,0), sendcount, MPI_DOUBLE, myRank-1, 4, &Hz(0,0,0), recvcount, MPI_DOUBLE, myRank-1, 2, MPI_COMM_WORLD, &status);
    }
}

void MPI_EH(int dummy){
    if (dummy & 0x1) MPI_E();
    if (dummy & 0x2) MPI_H();
}

void updateEH()
{



    updateH_main();
    updateH_PML();

    MPI_EH(2);

    updateE_main();
    updateE_PML();

    MPI_EH(1);


}
