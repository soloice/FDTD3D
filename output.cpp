#include "fdtd.h"
using namespace std;
void output_gather()
{
    int it;
    FILE *fp;
    char filename[30];
    char dir[MAX_NAME_LEN] = "./";
    sprintf(filename, "gather_%02d.dat", myRank);
    fp=fopen(strcat(dir, filename), "w+");
    //ofstream fout("./gather.dat");
    if(fp==NULL)
    {
        printf("File open error in output_gather!");
	exit(1);
    }
    for(ii=0; ii<nrec; ++ii)
    {
        for(it=0; it<nt; ++it)
        {
            fprintf(fp, "%e ", gather(ii,it));

        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void output_slice(int it,double *EH_total)
{

    char xfilename[30],yfilename[30],zfilename[30];
    char xdir[MAX_NAME_LEN] = "./Output/";
    char ydir[MAX_NAME_LEN] = "./Output/";
    char zdir[MAX_NAME_LEN] = "./Output/";
    
    sprintf(xfilename, "xSlice%05d_%d.dat", it,NUM_OF_PROCESS);
    sprintf(yfilename, "ySlice%05d_%d.dat", it,NUM_OF_PROCESS);
    sprintf(zfilename, "zSlice%05d_%d.dat", it,NUM_OF_PROCESS);
    //ofstream fout(filename);
    FILE *fp;
    fp=fopen(strcat(xdir, xfilename), "w+");
    //fp=fopen(xfilename, "w+");
    if(fp==NULL)
    {
        printf("File open error in output_slice!");
	exit(1);
    }
    if(nxslice != 0){
    for(i=0;i<nxslice;++i){
        for(j=0;j<ny;++j){
            for(k=0;k<nz;++k){
                fprintf(fp, "%e ", EH_total(slicex[i]+order,j,k));
            }}}}
    fclose(fp);

    fp=fopen(strcat(ydir, yfilename), "w+");
    if(fp==NULL)
    {
        printf("File open error");
    }
    if(nyslice !=0){
    for(j=0;j<nyslice;++j){
       for(i=0;i<nx;++i){
            for(k=0;k<nz;++k){
                fprintf(fp, "%e ", EH_total(i+order,slicey[j],k));
            }}}}
    fclose(fp);

    fp=fopen(strcat(zdir, zfilename), "w+");
    if(fp==NULL)
    {
        printf("File open error");
    }
    if(nzslice != 0){
    for(k=0;k<nzslice;++k){
        for(i=0;i<nx;++i){
            for(j=0;j<ny;++j){
                fprintf(fp, "%e ", EH_total(i+order,j,slicez[k]));
            }}}}
    fclose(fp);

}
void output_wavefield(int it)
{
    char filename[30];
    char dir[MAX_NAME_LEN] = "./Output/";
    sprintf(filename, "Wavefield%05d.dat", it);
    FILE *fp;
    fp=fopen(strcat(dir,filename), "w+");

    if(fp==NULL)
    {
        printf("File open error in output_wavefield!");
	exit(1);
    }
    for(i=0;i<nx;i+=output_step_x_of_wavefield)
    {
        for(j=0;j<ny;j+=output_step_x_of_wavefield)
        {
            for(k=0;k<nz;k+=output_step_x_of_wavefield)
            {
                fprintf(fp, "%e ", Ex(i,j,k));
            }
        }
    }
    fclose(fp);
}

