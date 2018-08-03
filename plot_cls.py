import numpy as np
import matplotlib.pyplot as plt
import in_out as io

nbins = 6
lmax = 5000

path_yes = 'outputs_may14_magnification_yes/run_fidcl.dat'
path_no = 'outputs_may14_magnification_no/run_fidcl.dat'

dic_yes = io.read_cls_class(path_yes)
dic_no = io.read_cls_class(path_no)

cls_dd_yes = dic_yes['cl_dd']
cls_dl_yes = dic_yes['cl_dl']
cls_ll_yes = dic_yes['cl_ll']
cls_dd_no = dic_no['cl_dd']
cls_dl_no = dic_no['cl_dl']
cls_ll_no = dic_no['cl_ll']

cls_dd_mag = abs(cls_dd_yes - cls_dd_no)
cls_dl_mag = abs(cls_dl_yes - cls_dl_no)
cls_ll_mag = abs(cls_ll_yes - cls_ll_no)

larr = np.linspace( 0, lmax, num = int(lmax+1) )

TITLE_SIZE = 30
LABEL_SIZE = 20
HANDLELENGTH = 4
COLUMNSPACING = 4
plt.rc('axes', titlesize=TITLE_SIZE)
plt.rc('axes', labelsize=TITLE_SIZE)



# dl
f, axarr = plt.subplots(nbins, nbins, sharex = True, sharey = True, figsize = (30,20))
for ii in range(0, nbins):
    for jj in range(0, nbins):
        not_mag = cls_dl_no[:,ii,jj] * larr**2
        mag = cls_dl_mag[:,ii,jj] * larr**2
        axarr[ii, jj].loglog(larr, not_mag, color = 'black', label = '$\ell^2C(\ell) (\mathsf{without\ magnification})$')
        axarr[ii, jj].loglog(larr, mag, '--', color = 'red', label = '$\mathsf{abs}\left[\ell^2C(\ell)\left(\mathsf{with\ magnification} - \mathsf{without\ magnification}  \\right)\\right]$')
        axarr[ii, jj].text(1, 0, '(' + str(ii+1) + ',' + str(jj+1) + ')', ha = 'right', va = 'bottom', transform=axarr[ii,jj].transAxes, fontsize = LABEL_SIZE)
        axarr[ii, jj].tick_params(which = 'both', direction = 'inout', bottom = True, top = True, left = True, right = True)
        if ii == nbins - 1 and jj == 0:
            axarr[ii, jj].legend(loc = 'center left', bbox_to_anchor=(0., -0.3), bbox_transform = axarr[ii,jj].transAxes, ncol = 2, fontsize = LABEL_SIZE, handlelength = HANDLELENGTH, columnspacing = COLUMNSPACING)
#f.suptitle(' ( clustering , shear ) ')
f.text(0.5, 0.92, '(clustering, shear)', transform = f.transFigure, ha = 'center', fontsize = TITLE_SIZE)
f.subplots_adjust(hspace=0)
f.subplots_adjust(wspace=0)
plt.savefig('clustering_shear_cls.pdf', bbox_inches='tight')

# dd
f, axarr = plt.subplots(nbins, nbins, sharex = True, sharey = True, figsize = (30,20))
for ii in range(0, nbins):
    for jj in range(0, nbins):
        not_mag = cls_dd_no[:,ii,jj] * larr**2
        mag = cls_dd_mag[:,ii,jj] * larr**2
        axarr[ii, jj].loglog(larr, not_mag, color = 'black', label = '$\ell^2C(\ell) (\mathsf{without\ magnification})$')
        axarr[ii, jj].loglog(larr, mag, '--', color = 'red', label = '$\mathsf{abs}\left[\ell^2C(\ell)\left(\mathsf{with\ magnification} - \mathsf{without\ magnification}  \\right)\\right]$')
        axarr[ii, jj].text(1, 0, '(' + str(ii+1) + ',' + str(jj+1) + ')', ha = 'right', va = 'bottom', transform=axarr[ii,jj].transAxes, fontsize = LABEL_SIZE)
        axarr[ii, jj].tick_params(which = 'both', direction = 'inout', bottom = True, top = True, left = True, right = True)
        if ii == nbins - 1 and jj == 0:
            axarr[ii, jj].legend(loc = 'center left', bbox_to_anchor=(0., -0.3), bbox_transform = axarr[ii,jj].transAxes, ncol = 2, fontsize = LABEL_SIZE, handlelength = HANDLELENGTH, columnspacing = COLUMNSPACING)
#f.suptitle(' ( clustering , shear ) ')
f.text(0.5, 0.92, '(clustering, clustering)', transform = f.transFigure, ha = 'center', fontsize = TITLE_SIZE)
f.subplots_adjust(hspace=0)
f.subplots_adjust(wspace=0)
plt.savefig('clustering_clustering_cls.pdf', bbox_inches='tight')

