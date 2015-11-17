import numpy as np
import os as os
import sys as sys
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import copy

FS=16

class ParamFisher:
    """ Fisher matrix parameter """
    val=0.0
    dval=0.0
    prior=0.0
    name="str"
    label="$x$"
    isfree=False
    do_plot=True

    def __init__(self,val,dval,prior,name,label,isfree,do_plot):
        self.val=val
        self.dval=dval
        self.prior=prior
        self.name=name
        self.label=label
        self.isfree=isfree
        self.do_plot=do_plot

def find_param(param_list,name):
    index=0
    for par in param_list:
        if par.name==name :
            return index
        index+=1
    sys.exit("No parameter "+name)

def plot_fisher_single(params,name,fishermat,ax,fc,lw,ls,lc,fact_axis) :
    nb=128

    sigma_max=0
    for i in np.arange(len(fishermat)) :
        i1=find_param(params,name)
        covar_full=np.linalg.inv(fishermat[i])
        sigma=np.sqrt(covar_full[i1,i1])
        if sigma>=sigma_max :
            sigma_max=sigma
        x_arr=params[i1].val-4*sigma+8*sigma*np.arange(nb)/(nb-1.)
        p_arr=np.exp(-(x_arr-params[i1].val)**2/(2*sigma**2))
        ax.set_title("$\\sigma($"+params[i1].label+"$)=%.3lf$"%sigma)
        ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i])
    ax.set_xlim([params[i1].val-fact_axis*sigma_max,params[i1].val+fact_axis*sigma_max])
    ax.set_xlabel(params[i1].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(FS-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(FS-2)

def plot_fisher_two(params,name1,name2,fishermat,ax,fc,lw,ls,lc,fact_axis) :
    sig0_max=0
    sig1_max=0
    for i in np.arange(len(fishermat)) :
        i1=find_param(params,name1)
        i2=find_param(params,name2)
        covar_full=np.linalg.inv(fishermat[i])
        covar=np.zeros([2,2])
        covar[0,0]=covar_full[i1,i1]
        covar[0,1]=covar_full[i1,i2]
        covar[1,0]=covar_full[i2,i1]
        covar[1,1]=covar_full[i2,i2]
        sig0=np.sqrt(covar[0,0])
        sig1=np.sqrt(covar[1,1])

        if sig0>=sig0_max :
            sig0_max=sig0
        if sig1>=sig1_max :
            sig1_max=sig1

        w,v=np.linalg.eigh(covar)
        angle=180*np.arctan2(v[1,0],v[0,0])/np.pi
        a_1s=np.sqrt(2.3*w[0])
        b_1s=np.sqrt(2.3*w[1])
        a_2s=np.sqrt(6.17*w[0])
        b_2s=np.sqrt(6.17*w[1])

        centre=np.array([params[i1].val,params[i2].val])

        e_1s=Ellipse(xy=centre,width=2*a_1s,height=2*b_1s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i],linestyle="solid",edgecolor=lc[i])
#                     facecolor=fc[i],linewidth=lw[i],linestyle=ls[i],edgecolor=lc[i])
        e_2s=Ellipse(xy=centre,width=2*a_2s,height=2*b_2s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i],linestyle="dashed",edgecolor=lc[i])
#                     facecolor=fc[i],linewidth=lw[i],linestyle=ls[i],edgecolor=lc[i])

        ax.add_artist(e_2s)
        ax.add_artist(e_1s)
        ax.set_xlim([params[i1].val-fact_axis*sig0_max,
                     params[i1].val+fact_axis*sig0_max])
        ax.set_ylim([params[i2].val-fact_axis*sig1_max,
                     params[i2].val+fact_axis*sig1_max])
        ax.set_xlabel(params[i1].label,fontsize=FS)
        ax.set_ylabel(params[i2].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(FS-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(FS-2)

def plot_fisher_all(params, #Parameters in the FMs
                    fishermat, #FMs to plot
                    fc,lw,ls,lc, #Foreground colors, line widths, line styles and line colours for each FM
                    labels, #Labels for each FM
                    fact_axis, #The x and y axes will be fact_axis x error in each parameter
                    fname) : #File to save the plot
    index_plot=np.where(np.array([p.do_plot for p in params]))
    n_params=len(index_plot[0])
    param_plot=params[index_plot]

    fig=plt.figure(figsize=(10,9))
    plt.subplots_adjust(hspace=0,wspace=0)
    for i in np.arange(n_params) : #Plot pdfs and ellipses
        i_col=i
        for j in np.arange(n_params-i)+i :
            i_row=j
            iplot=i_col+n_params*i_row+1

            ax=fig.add_subplot(n_params,n_params,iplot)
            if i==j :
                plot_fisher_single(params,param_plot[i].name,fishermat,
                                   ax,fc,lw,ls,lc,fact_axis)
            else :
                plot_fisher_two(params,param_plot[i].name,param_plot[j].name,
                                fishermat,ax,fc,lw,ls,lc,fact_axis)

            if i_row!=n_params-1 :
                ax.get_xaxis().set_visible(False)

            if i_col!=0 :
                ax.get_yaxis().set_visible(False)

            if i_col==0 and i_row==0 :
                ax.get_yaxis().set_visible(False)
                
            ax.locator_params(nbins=6)

    if n_params>1 : #Add labels in a separate plot
        ax=fig.add_subplot(n_params,n_params,2)
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])
        for i in np.arange(len(labels)) :
            ax.plot([-1,1],[-3,-3],color=lc[i],linestyle=ls[i],
                    linewidth=lw[i],label=labels[i])
        ax.legend(loc='upper left',frameon=False,fontsize=FS)
        ax.axis('off')

    if fname!="none" :
        plt.savefig(fname,bbox_inches='tight')

    plt.show()

