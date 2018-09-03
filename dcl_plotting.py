import numpy as np
import os as os
import matplotlib.pyplot as plt
import csv
import in_out as io


length = 5001
larr = np.linspace(1, length, num = length)

logarr = np.logspace(0, np.log10(length), num = 50, dtype = 'int', base = 10)

test_path = 'cut_dcl/'

#dndn_cases = ['A', 'B']
#edn_cases = ['a', 'b', 'c', 'd']
#nrows_edn = [0,1,2,3,4,5]

dndn_cases = ['DA', 'DB']
edn_cases = ['a', 'b']
nrows_edn = [1, 2]

twonbins = 4

path_yes = 'outputs_aug19_magnification_yes/run_fidcl.dat'
path_no = 'outputs_aug19_magnification_no/run_fidcl.dat'

dic_yes = io.read_cls_class(path_yes)
dic_no = io.read_cls_class(path_no)

cls_dd_yes = dic_yes['cl_dd']
cls_dl_yes = dic_yes['cl_dl']
cls_ll_yes = dic_yes['cl_ll']
cls_dd_no = dic_no['cl_dd']
cls_dl_no = dic_no['cl_dl']
cls_ll_no = dic_no['cl_ll']

for ii in range(0, len(dndn_cases)):
    for jj in range(0, len(edn_cases)):
        if True:
        #if edn_cases[jj] == 'd':
            for kk in range(0, len(nrows_edn)):
                print dndn_cases[ii], edn_cases[jj], str(nrows_edn[kk])
                folder = test_path + dndn_cases[ii] + '_' + edn_cases[jj] + str(nrows_edn[kk])

                dcl = np.zeros([twonbins, twonbins, length])
                for ll in range(0, twonbins):
                    for mm in range(0, twonbins):
                        with open(folder + '/' + str(ll) + '_' + str(mm) + '.ell', 'rb') as f:
                            reader = csv.reader(f)
                            nn = 0
                            for row in reader:
                                dcl[ll, mm, nn] = float(row[-1])
                                nn += 1
                                if nn >= length:
                                    break
                
                f, axarr = plt.subplots(twonbins, twonbins, sharex = True, sharey = True, figsize = (30, 20))
                for xx in range(0, twonbins):
                    for yy in range(0, twonbins):
                        signal = (dcl[xx, yy, :] * larr**2)
                        where_cut = np.zeros(len(signal))
                        where_cut[np.where(abs(signal)>1e-15)] = 1
                        not_mag = np.zeros(len(larr))
                        if xx < twonbins/2 and yy < twonbins/2: # clustering - clustering
                            not_mag = cls_dd_no[:,xx,yy] * larr**2
                        if xx < twonbins/2 and yy > (twonbins/2-1): # clustering - shear
                            not_mag = cls_dl_no[:,xx,yy-twonbins/2] * larr**2
                        if xx > (twonbins/2-1) and yy > (twonbins/2-1): # shear - shear
                            not_mag = cls_ll_no[:,xx-twonbins/2,yy-twonbins/2] * larr**2
                        if xx > (twonbins/2-1) and yy < twonbins/2: # shear - clustering
                            not_mag = cls_dl_no[:,yy,xx-twonbins/2] * larr**2
                        if xx > (twonbins/2-1) and yy > (twonbins/2-1): # shear - shear
                            not_mag = cls_ll_no[:,xx-twonbins/2,yy-twonbins/2] * larr**2
                        #not_mag = not_mag * where_cut
                        axarr[xx, yy].loglog(larr, not_mag, color = 'blue')
                        if max(abs(signal)) > 1e-15:
                            axarr[xx, yy].loglog(larr[np.where(signal < 0)], abs(signal[np.where(signal < 0)]), color = 'red', linestyle = 'dotted')
                            axarr[xx, yy].loglog(larr[np.where(signal > 0)], signal[np.where(signal > 0)], color = 'red', linestyle = 'solid')
                        axarr[xx, yy].set_ylim(bottom = 1e-8, top = 1e-1)
                        #axarr[xx, yy].set_xscale('symlog')
                        #axarr[xx, yy].set_yscale('symlog')
                f.subplots_adjust(hspace = 0)
                f.subplots_adjust(wspace = 0)
                f.text(0.5, 0.92, dndn_cases[ii] + '_' + edn_cases[jj] + '_' + str(nrows_edn[kk]), transform = f.transFigure, ha = 'center', fontsize = 30)
                plt.savefig(folder + '/' + dndn_cases[ii] + '_' + edn_cases[jj] + '_' + str(nrows_edn[kk]) + 'loglog.pdf', bbox_inches = 'tight')


        else:
            for kk in range(0, 1):
                folder = test_path + dndn_cases[ii] + '_' + edn_cases[jj] + str(nrows_edn[kk])

                dcl = np.zeros([twonbins, twonbins, length])
                for ll in range(0, twonbins):
                    for mm in range(0, twonbins):
                        with open(folder + '/' + str(ll) + '_' + str(mm) + '.ell', 'rb') as f:
                            reader = csv.reader(f)
                            nn = 0
                            for row in reader:
                                dcl[ll, mm, nn] = float(row[-1])
                                nn += 1
                                if nn >= length:
                                    break
                
                f, axarr = plt.subplots(twonbins, twonbins, sharex = True, sharey = True, figsize = (30, 20))
                ####
                for xx in range(0, twonbins):
                    for yy in range(0, twonbins):
                        signal = (dcl[xx, yy, :] * larr**2)
                        where_cut = np.zeros(len(signal))
                        where_cut[np.where(abs(signal)>1e-15)] = 1
                        not_mag = np.zeros(len(larr))
                        if xx < twonbins/2 and yy < twonbins/2: # clustering - clustering
                            not_mag = cls_dd_no[:,xx,yy] * larr**2
                        if xx < twonbins/2 and yy > (twonbins/2-1): # clustering - shear
                            not_mag = cls_dl_no[:,xx,yy-twonbins/2] * larr**2
                        if xx > (twonbins/2-1) and yy > (twonbins/2-1): # shear - shear
                            not_mag = cls_ll_no[:,xx-twonbins/2,yy-twonbins/2] * larr**2
                        if xx > (twonbins/2-1) and yy < twonbins/2: # shear - clustering
                            not_mag = cls_dl_no[:,yy,xx-twonbins/2] * larr**2
                        if xx > (twonbins/2-1) and yy > (twonbins/2-1): # shear - shear
                            not_mag = cls_ll_no[:,xx-twonbins/2,yy-twonbins/2] * larr**2
                        #not_mag = not_mag * where_cut
                        axarr[xx, yy].loglog(larr, not_mag, color = 'blue')
                        if max(abs(signal)) > 1e-15:
                            axarr[xx, yy].loglog(larr[np.where(signal < 0)], abs(signal[np.where(signal < 0)]), color = 'red', linestyle = 'dotted')
                            axarr[xx, yy].loglog(larr[np.where(signal > 0)], signal[np.where(signal > 0)], color = 'red', linestyle = 'solid')
                        axarr[xx, yy].set_ylim(bottom = 1e-8, top = 1e-1)
                        #axarr[xx, yy].set_xscale('symlog')
                        #axarr[xx, yy].set_yscale('symlog')
                ####
                f.subplots_adjust(hspace = 0)
                f.subplots_adjust(wspace = 0)
                f.text(0.5, 0.92, dndn_cases[ii] + '_' + edn_cases[jj] + '_' + str(nrows_edn[kk]), transform = f.transFigure, ha = 'center', fontsize = 30)
                plt.savefig(folder + '/' + dndn_cases[ii] + '_' + edn_cases[jj] + '_' + str(nrows_edn[kk]) + 'loglog.pdf', bbox_inches = 'tight')
                #for xx in range(0, twonbins):
                #    for yy in range(0, twonbins):
                #        signal = abs(dcl[xx, yy, :] * larr**2)
                #        if max(signal) < 1e-15:
                #            continue
                #        signal[np.where(signal < 1e-15)] = 0.5*min(signal[np.where(signal > 1e-15)])
                #        axarr[xx, yy].loglog(larr[logarr], signal[logarr])
                #f.subplots_adjust(hspace = 0)
                #f.subplots_adjust(wspace = 0)
                #f.text(0.5, 0.92, dndn_cases[ii] + '_' + edn_cases[jj] + '_' + str(nrows_edn[kk]), transform = f.transFigure, ha = 'center', fontsize = 30)
                #plt.savefig(folder + '/' + dndn_cases[ii] + '_' + edn_cases[jj] + '_' + str(nrows_edn[kk]) + '.pdf', bbox_inches = 'tight')
