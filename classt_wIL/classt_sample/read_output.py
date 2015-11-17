import numpy as np
import struct as struct
import sys as sys
import matplotlib.pyplot as plt

def read_cls_dd_binary(fname) :
    f=open(fname,"rd")
    data=f.read()
    f.close()
    
    lmax,nbins=struct.unpack('2I',data[0:8])

    cl_out=np.zeros([lmax+1,nbins,nbins])
    ncross=0
    sizecross=8+(lmax-1)*8
    for i in np.arange(nbins) :
        for j in np.arange(nbins-i)+i :
            offset=8+ncross*sizecross
            i1,i2=struct.unpack('2I',data[offset:offset+8])
            if (i1!=i) and (i2!=j) :
                print "Error reading file"
                sys.exit(1)
            cl=struct.unpack('%dd'%(lmax-1),data[offset+8:offset+sizecross])
            cl_out[2:,i,j]=cl
            if i!=j :
                cl_out[:,j,i]=cl_out[:,i,j]
            ncross+=1
            
    return lmax,nbins,cl_out

def read_cls_dd_ascii(fname) :
    data=np.loadtxt(fname,unpack=True)
    ncross=len(data)-1
    lmax=len(data[0])+1
    nbins=int((-1+np.sqrt(1+8*ncross))/2)
    if nbins*(nbins+1)/2!=ncross :
        print "Error reading file"
        sys.exit(1)
    cl_out=np.zeros([lmax+1,nbins,nbins])
    ncross=0
    for i in np.arange(nbins) :
        for j in np.arange(nbins-i)+i :
            cl_out[2:,i,j]=data[ncross+1,:]
            if i!=j :
                cl_out[:,j,i]=cl_out[:,i,j]
            ncross+=1
            
    return lmax,nbins,cl_out

def read_cls_dd(fname,use_binary) :
    if use_binary==1 :
        return read_cls_dd_binary(fname)
    else :
        return read_cls_dd_ascii(fname)

lmax,nbins,cl1=read_cls_dd('out_cl_dd.dat',1)
lmax,nbins,cl2=read_cls_dd('out_cl_dd.txt',0)
larr=np.arange(lmax-1)+2
for i in np.arange(nbins) :
    for j in np.arange(nbins-i)+i :
        plt.title("%d "%(i+1)+"%d"%(j+1))
        plt.plot(larr,np.fabs(cl1[2:,i,j]-cl2[2:,i,j]))
#        plt.plot(larr,cl1[2:,i,j]*2*np.pi/(larr*(larr+1)),'r',linewidth=3)
#        plt.plot(larr,cl2[2:,i,j]*2*np.pi/(larr*(larr+1)),'b')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.show()
