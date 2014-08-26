#!/usr/bin/python
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pylab as pylab
import os
import time
import re

def read_slice(fname):
    with open(fname) as fslice:
        slice_nx,slice_ny,slice_nz = fslice.readline().split()
        slice_x = fslice.readline().split()
        slice_y = fslice.readline().split()
        slice_z = fslice.readline().split()
        slice_nx = int(slice_nx);slice_ny = int(slice_ny);slice_nz = int(slice_nz)
        return slice_nx,slice_ny,slice_nz

def read_rec(frec):
    global nrec
    with open(frec) as fp:
        nrec = int(file.readline().strip('\n'))



def read_par():
    global nx,ny,nz,slice_nx,slice_ny,slice_nz,nt
    with open('par.in') as fpar:
        fpar.readline()
        dx,dy,dz,dt = fpar.readline().split()
        print 'dx dy dz dt: ',dx,dy,dz,dt
        fpar.readline()
        nx,ny,nz,nt = fpar.readline().split()
        nx = int(nx);ny = int(ny);nz = int(nz);nt=int(nt)
        print 'nx ny nz nt: ',nx,ny,nz,nt
        fpar.readline()
        nt_src = fpar.readline()
        print 'nt of src: ',nt_src
        fpar.readline()
        step_t_wavefield,step_x_wavefield = fpar.readline().split()
        print 'output time step and space step of wavefidld: ',step_t_wavefield,step_x_wavefield
        fpar.readline()
        step_slice = fpar.readline()
        print 'output step of slice: ',step_slice
        fpar.readline()
        npml_x,npml_y,npml_z= fpar.readline().split()
        print 'npml x y z: ',npml_x,npml_y,npml_z
        fpar.readline()
        fpar.readline() #pml m kapxmax kapymax kapzmax alpha
        fpar.readline()
        fsrc= fpar.readline().strip('\n')
        print 'src.in: ',fsrc
        fpar.readline()
        frec= fpar.readline().strip('\n')
        print 'rec.in: ',frec
        fpar.readline()
        feps = fpar.readline().strip('\n')
        fpar.readline()
        fmu = fpar.readline().strip('\n')
        fpar.readline()
        fsig= fpar.readline().strip('\n')
        fpar.readline()
        fslice= fpar.readline().strip('\n')
        slice_nx,slice_ny,slice_nz = read_slice(fslice)





def view_slice():
    xlist = os.popen('ls xSlice*dat').readlines()
    ylist = os.popen('ls ySlice*dat').readlines()
    zlist = os.popen('ls zSlice*dat').readlines()
    i = 0
    for xname in xlist:
        print xname
        yname = ylist[i]
        zname = zlist[i]
        i += 1
        xdata = loadtxt(xname.strip('\n'))
        ydata = loadtxt(yname.strip('\n'))
        zdata = loadtxt(zname.strip('\n'))

        xslice = reshape(xdata,(slice_nx,ny,nz))
        yslice = reshape(ydata,(slice_ny,nx,nz))
        zslice = reshape(zdata,(slice_nz,nx,ny))
        # data = reshape(data,(126,101))
        clf()
        imshow(xslice[0])
        colorbar()
        savefig(re.findall("^\w+",xname)[0]+".jpg")
        clf()
        imshow(yslice[0])
        colorbar()
        savefig(re.findall("^\w+",yname)[0]+".jpg")
        clf()
        imshow(zslice[0])
        colorbar()
        savefig(re.findall("^\w+",zname)[0]+".jpg")


#        show()
#        show(block=False)
#        time.sleep(0.5)
#        close()

def view_gather():
    global  nrec,nt
    ilist = os.popen('ls gather*dat').readlines()
    i = 0
    for name in ilist:
        gather = loadtxt(name.strip('\n'))
        if gather.max() == 0 and gather.min() == 0:
            continue
    # for i in range(len(gather)):
    	plot(gather[:]/max(abs(gather[:]))+i)
        i += 1
    savefig('gather.png')


def check():

    data1 = loadtxt("ySlice00014_2.dat")
    # data2 = loadtxt("ySlice00014_3.dat")


    slice1 = reshape(data1,(slice_ny,nx,nz))
    # slice2 = reshape(data2,(slice_ny,nx,nz))
    # slice_diff = slice1 - slice2;
    # zslice = reshape(zdata,(slice_nz,nx,ny))
    # data = reshape(data,(126,101))
    figure()
    imshow(slice1[0])
    colorbar()
    # # savefig(re.findall("^\w+",xname)[0]+".jpg")
    # figure()
    # imshow(slice2[0])
    # colorbar()
    # # savefig(re.findall("^\w+",yname)[0]+".jpg")
    # figure()
    # imshow(slice_diff[0])
    # colorbar()
    show()
    # imshow()
    # savefig(re.findall("^\w+",zname)[0]+".jpg")

# read_par()
# os.chdir("./Output/")
# check()


        







read_par()
view_gather()
os.chdir("./Output/")
view_slice()

#view_wavefield()

    
