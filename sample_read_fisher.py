import numpy as np
import matplotlib.pyplot as plt

def read_fisher(prefix) :
    data=np.load(prefix+"/fisher_raw.npz")
    
    fish_package={'names' :data['names'],
                  'values':data['values'],
                  'labels':data['labels'],
                  'fisher':data['fisher_tot']}
    return fish_package

def print_sigmas(prefix) :
    pkg=read_fisher(prefix)
    cov=np.linalg.inv(pkg['fisher'])
    for i,nm in enumerate(pkg['names']) :
        print nm+" = %lE"%(pkg['values'][i])+" +- %lE"%(np.sqrt(cov[i,i]))

print_sigmas("sample_files")

