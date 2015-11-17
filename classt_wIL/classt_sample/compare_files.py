import struct
import numpy as np
import matplotlib.pyplot as plt
import time

nopen=5000

print "Reading ascii"
start=time.time()
for i in np.arange(nopen) :
    larr_1,dlarr_1=np.loadtxt("double_old_cl_d1d1.dat",unpack=True)
end=time.time()
print end-start

print "Reading binary"
start=time.time()
for i in np.arange(nopen) :
    f=open("double_old_b_cl_d1d1.dat")
    f.seek(4)
    dlarr_2=np.fromfile(f,dtype='d64')
end=time.time()
print end-start

print "done"
'''
for b1 in np.arange(bmax)+1 :
    for b2 in np.arange(bmax-b1+1)+b1 :
        print b1,b2
        fname_1="double_old_cl_d%d"%b1+"d%d.dat"%b2
        fname_2="double_old_b_cl_d%d"%b1+"d%d.dat"%b2

        larr_1,dlarr_1=np.loadtxt(fname_1,unpack=True)

        f=open(fname_2)
        lmax=struct.unpack('i',f.read(4))[0]
        larr_2=np.arange(lmax-1)+2
        dlarr_2=np.fromfile(f,dtype='d64')
        
        plt.plot(larr_1,(dlarr_1-dlarr_2)/larr_1/(larr_1+1))
        plt.show()
'''
