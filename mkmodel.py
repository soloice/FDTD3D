#!/usr/bin/python
from numpy import *
from scipy import signal
from matplotlib.pyplot import *
import os

def finddt(epmin,mumin,dx,dy,dz):

    epmin = epmin*ep0;
    mumin = mumin*mu0;
    dtmax = 6.0/7.0*sqrt(epmin*mumin/(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2));
    return dtmax
    
def finddx(epmax,mumax,fmax):
    mu0 = 1.2566370614 * 10**-6;
    ep0 = 8.8541878176 * 10**-12;
    epmax = epmax*ep0;
    mumax = mumax*mu0;
    wlmin = 1/(fmax*sqrt(epmax*mumax));
    dxmax = wlmin/8.0;
    return dxmax;

def findnt(epmax,mumax,dt,Z):
    mu0 = 1.2566370614 * 10**-6;
    ep0 = 8.8541878176 * 10**-12;
    epmax = epmax*ep0;
    mumax = mumax*mu0;
    nt = int(1.0/2.0 * ceil(Z * sqrt(epmax*mumax) / dt));
    return nt
    
def blackharrispulse(fmax,dt,nt):
    a = [0.35322222,-0.488,0.145,-0.010222222]
    T = 1.14/fmax
    # t = arange(0,T,dt)
    t = arange(nt+1)*dt
    window = zeros(size(t));
    # window = zeros(1251)
    for n in range(4):
        window = window + a[n]*cos(2*n*pi*t/T)

    window[t>=T] = 0
    p = window;
    p = append(window[2:],0) - window[1:];
    p = p/max(abs(p))
    return p

#######################################################



#%%
def discretization():
    global dx,dy,dz,dt,nx,ny,nz,nt,nt_src,srcpulse
    # dx = finddx(epmax,mumax,fmax);
    dx = 0.04;
    dy = dx; dz = dx;

    print "dx dy dz: ",dx,dy,dz;
    
    # nx = int(ceil(X/dx)); ny = int(ceil(Y/dy)); nz = int(ceil(Z/dz));
    nx = 99
    ny = 85
    nz = 70
    print "nx ny nz: ",nx,ny,nz;

    
    # dt = finddt(epmin,mumin,dx,dy,dz)
    dt = 8.0e-11
    print "dt: ",dt

    
    # nt = findnt(epmax,mumax,dt,Z)
    nt = 651

    srcpulse = blackharrispulse(fmax,dt,nt)


    nt_src = len(srcpulse)
    print "nt, nt_src: ",nt,nt_src

def src_rec():
    global src,rec
    #nxprop = nx * 2 -1;nyprop = ny * 2 -1;nzprop = nz * 2 -1;
    # nsrcx = 1;nsrcy = 1;nsrcz = 1
    nrecx = 1;nrecy = 1;nrecz = 1
    # src = zeros((nsrcx*nsrcy*nsrcz,3)); rec = zeros((nrecx*nrecy*nrecz,3))
    nsrc = 8
    src = []
    src.append([24,30,30,'Hy',1]);src.append([24,31,30,'Hy',1]);#src.append([72,73,74,'Hz',1]);
    src.append([27,30,30,'Hy',-1]);src.append([27,31,30,'Hy',-1]);#src.append([74,73,74,'Hz',-1]);
    src.append([25,32,30,'Hx',1]);src.append([26,32,30,'Hx',1]);#src.append([74,73,74,'Hx',1]);
    src.append([25,29,30,'Hx',-1]);src.append([26,29,32,'Hx',-1]);#src.append([74,73,72,'Hx',-1]);
    rec = zeros((nrecx*nrecy*nrecz,3));rec[0,:] = 30

    # for i in range(nsrcx):
    #     for j in range(nsrcy):
    #         for k in range(nsrcz):
    #             src[k + nsrcz*j + nsrcy*nsrcz*i,0] = fix(((i+1) * X/(nsrcx+1.0))/dx);
    #             src[k + nsrcz*j + nsrcy*nsrcz*i,1] = fix(((j+1) * Y/(nsrcy+1.0))/dy);
    #             src[k + nsrcz*j + nsrcy*nsrcz*i,2] = int(Zsurf/dx);
                
    # for i in range(nrecx):
    #     for j in range(nrecy):
    #         for k in range(nrecz):
    #             rec[k + nrecz*j + nrecy*nrecz*i,0] = fix(((i+1) * X/(nrecx+1.0))/dx);
    #             rec[k + nrecz*j + nrecy*nrecz*i,1] = fix(((j+1) * Y/(nrecy+1.0))/dy);
    #             rec[k + nrecz*j + nrecy*nrecz*i,2] = int(Zsurf/dx);
    print "nsrc: ",nsrc

    
    #%%write data
    with file('src.in', 'w') as fsrc:
        # fsrc.write("%d %d\n" %(nsrcx*nsrcy*nsrcz,nt_src))
#        savetxt(fsrc,src,fmt='%d')
        fsrc.write("%d %d\n" %(nsrc,nt_src))
        for i in range(nsrc):
            fsrc.write("%d %d %d %s\n" %(src[i][0],src[i][1],src[i][2],src[i][3]))
        for i in range(nsrc):
            savetxt(fsrc,srcpulse*src[i][4])
    
    print "nrec: ",nrecx*nrecy*nrecz
    with file('rec.in', 'w') as frec:
        frec.write("%d\n" %(nrecx*nrecy*nrecz))
#        savetxt(frec,rec,fmt='%d')
        for i in range(nrecx*nrecy*nrecz):
            frec.write("%d %d %d %s\n" %(rec[i,0],rec[i,1],rec[i,2],"Ez"))

def islice():
    nslicex = 1;nslicey = 1;nslicez = 1;
    # slicex = fix(linspace(0,X,nslicex+1,endpoint=False)[1:]/dx)
    # slicey = fix(linspace(0,X,nslicey+1,endpoint=False)[1:]/dy)
    # slicez = fix(linspace(0,X,nslicez+1,endpoint=False)[1:]/dz)
    slicex = [];slicey = [];slicez = [];
    slicex.append([24,'Ez'])
    slicey.append([30,'Ez'])
    slicez.append([30,'Ez'])
    print "nslicex,nslicey,nslicez: ",nslicex,nslicey,nslicez
    
    with file('slice.in', 'w') as fslice:
        fslice.write("%d %d %d\n" % (len(slicex),len(slicey),len(slicez)))
        for i in range(len(slicex)):
            fslice.write("%d %s\n" %(slicex[0][0],slicex[0][1]))
        for i in range(len(slicey)):
            fslice.write("%d %s\n" %(slicey[0][0],slicey[0][1]))
        for i in range(len(slicez)):
            fslice.write("%d %s\n" %(slicez[0][0],slicez[0][1]))

def distance(x,y,z,x0,y0,z0):
    return sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)    
        
def ep_mu_sig():
    global ep_plot,mu_plot,sig_plot
    ep_plot=zeros((ny,nz))
    mu_plot=zeros((ny,nz))
    sig_plot=zeros((ny,nz))
    x0 = 2.5; y0 = 2.5; z0 = 2.5;
    
    with file('eps.in', 'w') as fep:
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if i >= 20 and i <=30 and j >= 25 and j <= 35 and k == 90:
                        fep.write("%f " % (15.0))
                    else:
                        fep.write("%f " % (9.0))

                    # if i == 23: ep_plot[j,k] = 9.0
                    # idist = distance(i*dx,j*dy,k*dz,x0,y0,z0)
#                    if k <= int(Zsurf/dz):
#                        fep.write("%f " % (1.0))
#                        if i == int(nx/2): ep_plot[j,k] = 1.0
#                    elif idist < 0.222:
#                    if idist < 0.222:
#                        fep.write("%f " % (12.0))
#                        if i == int(nx/2): ep_plot[j,k] = 12.0
#                    else:
#                        fep.write("%f " % (9.0))
#                        if i == int(nx/2): ep_plot[j,k] = 9.0
    #                    fep.write("%f " % (1.0))
    #                    ep_plot[j,z] = 1.0
    #                elif z > int(nz/2.0) and j <= int(ny/2.0):
    #                    fep.write("%f " % (9.0))
    #                    ep_plot[j,z] = 9.0
    #                elif z > int(nz/3.0) and j > int(ny/2.0):
    #                    fep.write("%f " % (10.0))
    #                    ep_plot[j,z] = 10.0
    #                else:
    #                    fep.write("%f " % (8.0))
    #                    ep_plot[j,z] = 8.0
    #
    #                 if k < int(nz/2.0):
    #                     fep.write("%f " % (6.0))
    #                     if i == int(nx/2): ep_plot[j,k] = 6.0
    #                 else:
    #                     fep.write("%f " % (8.0))
    #                     if i == int(nx/2): ep_plot[j,k] = 8.0

               
                    
    with file('mu.in', 'w') as fmu:
        for i in range(nx):
            for j in range(ny):
                for z in range(nz):
                    fmu.write("%f " %(1.0))
                    mu_plot[j,z]=1.0
                    
    with file('sig.in', 'w') as fsig:
        for i in range(nx):
            for j in range(ny):
                for z in range(nz):
                    fsig.write("%e " %(1.0e-3))
                    sig_plot[j,z]=0
                

#%%
def par():
    with file('par.in', 'w') as fpar:
        fpar.write("#dx dy dz dt\n")
        fpar.write("%e %e %e %e\n" %(dx,dy,dz,dt))
        fpar.write("#nx ny nz nt\n")
        fpar.write("%d %d %d %d\n" %(nx,ny,nz,nt))
        fpar.write("#nt of src\n")
        fpar.write("%d\n" %(nt_src))
        fpar.write("#output time step and space step of wavefield\n")
        fpar.write("%d %d\n" % (outstep_t_wavefield, outstep_x_wavefield))
        fpar.write("#output step of slice\n")
        fpar.write("%d\n" % (outstep_slice))
        fpar.write("#npml x y z\n")
        fpar.write("12 12 12\n")
        fpar.write("#pml m kapxmax kapymax kapzmax alpha\n")
        fpar.write("4 5 5 5 0.0\n")
        fpar.write("#location of src\n")
        fpar.write("src.in\n")
        fpar.write("#location of rec\n")
        fpar.write("rec.in\n")
        fpar.write("#epsilon file\n")
        fpar.write("eps.in\n")
        fpar.write("#mu file\n")
        fpar.write("mu.in\n")
        fpar.write("#sigma file\n")
        fpar.write("sig.in\n")
        fpar.write("#slices file\n")
        fpar.write("slice.in\n")


#%%plot data

def iplot():
    """plot dat"""
    figure()
    plot(srcpulse)


    figure()
    subplot(1,3,1)
    imshow(ep_plot.T)
    
    subplot(1,3,2)
    imshow(mu_plot.T)
    
    subplot(1,3,3)
    imshow(sig_plot.T)
    
    
    figure()
    subplot(2,1,1)
    hold(True)
    for i in range(len(src[:,0])):
        plot(src[i,0]*dx,src[i,1]*dy,'*')
    axis('equal')
    # xlim(0,X);ylim(0,Y)

    subplot(2,1,2)
    hold(True)
    for i in range(len(rec[:,0])):
        plot(rec[i,0]*dx,rec[i,1]*dy,'*')
    axis('equal')
    # xlim(0,X);ylim(0,Y)

    show()


mu0 = 1.2566370614 * 10 ** -6;
ep0 = 8.8541878176 * 10 ** -12;
X = 5; Y = 5; Z = 4; #m
Zsurf = 2;#m
epmin = 1.0; mumin = 1.0;
epmax = 81.0; mumax = 1.1;
fmax = 100e6; #Hz
outstep_t_wavefield = 100000
outstep_x_wavefield = 10 
outstep_slice = 20

discretization()
src_rec()
ep_mu_sig()
islice()
par()
idir = './Output/' 
if not os.path.exists(idir):
    os.mkdir(idir)     
# iplot()
        
    
"""
np.savetxt('ep.in',ep,fmt='%lf')
np.savetxt('mu.in',mu,fmt='%lf')
np.savetxt('sig.in',sig,fmt='%lf')
"""

