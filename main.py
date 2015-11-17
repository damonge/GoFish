import numpy as np
import matplotlib.pyplot as plt
import common as com
import sys as sys
import os as os

if len(sys.argv)!=2:
    sys.exit("Usage: python main.py fname_params")
fname_params=sys.argv[1]
#sys.stderr=open('errorlog.txt','w')
print " "
print "__________________________"
print "| ------- GoFish ------- |"
print "|                        |"
print "|   _///_       _///_    |" 
print "|  /o    \/    /o    \/  |"   
print "|  > ))_./\    > ))_./\  |"   
print "|     <           <      |"   
print "|------------------------|"
print " "
print " "

print "<> Reading parameters"
par=com.ParamRun(fname_params)
print " "

print "<> Computing/reading relevant signal power spectra"
par.get_cls_all()
print " "

print "<> Computing relevant noise power spectra"
par.get_cls_noise()
par.plot_cls()
print " "

print "<> Computing Fisher matrix"
par.get_fisher_l()
par.plot_fisher()
print " "

#sys.stderr.close()
#sys.stderr=sys.__stderr__
#print "<> Contents of the error log :"
#print "### "
#print " "
#os.system("cat errorlog.txt")
#os.system("rm errorlog.txt")
#print " "
#print "### "
#print " "

print " "
print "__________________________"
print "|          Done!         |"
print "|------------------------|"
print " "

