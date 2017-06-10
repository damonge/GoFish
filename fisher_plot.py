import numpy as np
import os as os
import sys as sys
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import copy
from decimal import Decimal

FS=16

class ParamFisher:
    """ Fisher matrix parameter """
    val=0.0
    dval=0.0
    onesided=0
    name="str"
    label="$x$"
    isfree=False
    do_plot=True

    def __init__(self,val,dval,name,label,isfree,do_plot,onesided):
        self.val=val
        self.dval=dval
        self.name=name
        self.label=label
        self.isfree=isfree
        self.do_plot=do_plot
        self.onesided=onesided

def find_param(param_list,name):
    index=0
    for par in param_list:
        if par.name==name :
            return index
        index+=1
    sys.exit("No parameter "+name)

def plot_fisher_single(params,name,fishermat,ax,fc,lw,ls,lc,fact_axis,show_title=True, legend=False, labels=[]) :
    nb=128

    sigma_max=0
    i1=find_param(params,name)
    for i in np.arange(len(fishermat)) :
        covar_full=np.linalg.inv(fishermat[i])
        sigma=np.sqrt(covar_full[i1,i1])
        if sigma>=sigma_max :
            sigma_max=sigma
        x_arr=params[i1].val-4*sigma+8*sigma*np.arange(nb)/(nb-1.)
        p_arr=np.exp(-(x_arr-params[i1].val)**2/(2*sigma**2))
        if sigma < 1e-2:
            sigma_str = "$"+r'{} \times 10^{{{}}}$'.format(*('%.1E' % Decimal(sigma)).split('E'))
        else:
            sigma_str= "%.2lf"%sigma
            sigma_str= "%.2lf"%sigma
        if i == 0:
            # p_title = "$\\sigma($"+str(params[i1].label)+sigma_str
            p_title = sigma_str
        else:
            # p_title += "\n$\\sigma($"+str(params[i1].label)+sigma_str
            p_title += "\n"+sigma_str
        if labels != []:
            ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i], label=labels[i])
        else:
            ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i])
    if show_title:
        ax.set_title(p_title)
    ax.set_xlim([params[i1].val-fact_axis*sigma_max,params[i1].val+fact_axis*sigma_max])
    ax.set_xlabel(params[i1].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(FS-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(FS-2)
    

def plot_fisher_two(params,name1,name2,fishermat,ax,fc,lw,ls,lc,fact_axis) :
    sig0_max=0
    sig1_max=0
    i1=find_param(params,name1)
    i2=find_param(params,name2)
    for i in np.arange(len(fishermat)) :
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
                     facecolor=fc[i],linewidth=lw[i]/2.,linestyle="dashed",edgecolor=lc[i])
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
                    fname,  #File to save the plot
                    param_list=None, figsize=(20,20), add_label=True, show_title=True, label_loc=4) :
    
    if param_list is None:
        index_plot=np.where(np.array([p.do_plot for p in params]))
    else:
        index_plot_in = [find_param(params,name) for name in param_list]
        index_plot_toplot = list(np.where(np.array([p.do_plot for p in params]))[0])
        index_plot = [list(set(index_plot_in).intersection(set(index_plot_toplot)))]
    
    n_params=len(index_plot[0])
    param_plot=params[index_plot]

    fig=plt.figure(figsize=figsize)
    plt.subplots_adjust(hspace=0,wspace=0)
    for i in np.arange(n_params) : #Plot pdfs and ellipses
        i_col=i
        for j in np.arange(n_params-i)+i :
            i_row=j
            iplot=i_col+n_params*i_row+1

            ax=fig.add_subplot(n_params,n_params,iplot)
            if i==j :
                plot_fisher_single(params,param_plot[i].name,fishermat,
                                   ax,fc,lw,ls,lc,fact_axis,show_title=show_title)
            else :
                plot_fisher_two(params,param_plot[i].name,param_plot[j].name,
                                fishermat,ax,fc,lw,ls,lc,fact_axis)

            if i_row!=n_params-1 :
                ax.get_xaxis().set_visible(False)

            if i_col!=0 :
                ax.get_yaxis().set_visible(False)

            if i_col==0 and i_row==0 :
                ax.get_yaxis().set_visible(False)
                
            ax.locator_params(nbins=3)

    if add_label:
        if n_params>1 : #Add labels in a separate plot
            ax=fig.add_subplot(n_params,n_params,label_loc)
            ax.set_xlim([-1,1])
            ax.set_ylim([-1,1])
            
            for i in range(len(fishermat)) :
                ax.plot([-1,1],[-3,-3],color=lc[i],linestyle=ls[i],
                        linewidth=lw[i],label=labels[i])
            ax.legend(loc='upper left',frameon=False,fontsize=FS)
            ax.axis('off')

    if fname!="none" :
        plt.savefig(fname,bbox_inches='tight')

    plt.show()

def plot_fisher_single_params(params,name,fishermat,ax,fc,lw,ls,lc,fact_axis,show_title=True, legend=False, labels=[]) :
    nb=128

    for i in np.arange(len(fishermat)) :
        sigma_max=0
        i1=find_param(params[i],name)
        covar_full=np.linalg.inv(fishermat[i])
        sigma=np.sqrt(covar_full[i1,i1])
        if sigma>=sigma_max :
            sigma_max=sigma
        x_arr=params[i][i1].val-4*sigma+8*sigma*np.arange(nb)/(nb-1.)
        p_arr=np.exp(-(x_arr-params[i][i1].val)**2/(2*sigma**2))
        if sigma < 1e-2:
            sigma_str = "$"+r'{} \times 10^{{{}}}$'.format(*('%.1E' % Decimal(sigma)).split('E'))
        else:
            sigma_str= "%.2lf"%sigma
            sigma_str= "%.2lf"%sigma
        if i == 0:
            # p_title = "$\\sigma($"+str(params[i1].label)+sigma_str
            p_title = sigma_str
        else:
            # p_title += "\n$\\sigma($"+str(params[i1].label)+sigma_str
            p_title += "\n"+sigma_str
        if labels != []:
            ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i], label=labels[i])
        else:
            ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i])
    if show_title:
        ax.set_title(p_title)
    ax.set_xlim([params[i][i1].val-fact_axis*sigma_max,params[i][i1].val+fact_axis*sigma_max])
    ax.set_xlabel(params[i][i1].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(FS-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(FS-2)
    

def plot_fisher_two_params(params,name1,name2,fishermat,ax,fc,lw,ls,lc,fact_axis) :
    sig0_max=0
    sig1_max=0

    for i in np.arange(len(fishermat)) :
        i1=find_param(params[i],name1)
        i2=find_param(params[i],name2)
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

        centre=np.array([params[i][i1].val,params[i][i2].val])

        e_1s=Ellipse(xy=centre,width=2*a_1s,height=2*b_1s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i],linestyle="solid",edgecolor=lc[i])
#                     facecolor=fc[i],linewidth=lw[i],linestyle=ls[i],edgecolor=lc[i])
        e_2s=Ellipse(xy=centre,width=2*a_2s,height=2*b_2s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i]/2.,linestyle="dashed",edgecolor=lc[i])
#                     facecolor=fc[i],linewidth=lw[i],linestyle=ls[i],edgecolor=lc[i])

        ax.add_artist(e_2s)
        ax.add_artist(e_1s)
        ax.set_xlim([params[i][i1].val-fact_axis*sig0_max,
                     params[i][i1].val+fact_axis*sig0_max])
        ax.set_ylim([params[i][i2].val-fact_axis*sig1_max,
                     params[i][i2].val+fact_axis*sig1_max])
        ax.set_xlabel(params[i][i1].label,fontsize=FS)
        ax.set_ylabel(params[i][i2].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(FS-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(FS-2)

def plot_fisher_all_params(params, #Parameters in the FMs
                    fishermat, #FMs to plot
                    fc,lw,ls,lc, #Foreground colors, line widths, line styles and line colours for each FM
                    labels, #Labels for each FM
                    fact_axis, #The x and y axes will be fact_axis x error in each parameter
                    fname,  #File to save the plot
                    param_list=None, figsize=(20,20), add_label=True, show_title=True, label_loc=4) :
        
    n_params=len(param_list)
    param_plot=param_list

    fig=plt.figure(figsize=figsize)
    plt.subplots_adjust(hspace=0,wspace=0)
    for i in np.arange(n_params) : #Plot pdfs and ellipses
        i_col=i
        for j in np.arange(n_params-i)+i :
            i_row=j
            iplot=i_col+n_params*i_row+1

            ax=fig.add_subplot(n_params,n_params,iplot)
            if i==j :
                plot_fisher_single_params(params,param_plot[i],fishermat,
                                   ax,fc,lw,ls,lc,fact_axis,show_title=show_title)
            else :
                plot_fisher_two_params(params,param_plot[i],param_plot[j],
                                fishermat,ax,fc,lw,ls,lc,fact_axis)

            if i_row!=n_params-1 :
                ax.get_xaxis().set_visible(False)

            if i_col!=0 :
                ax.get_yaxis().set_visible(False)

            if i_col==0 and i_row==0 :
                ax.get_yaxis().set_visible(False)
                
            ax.locator_params(nbins=3)

    if add_label:
        if n_params>1 : #Add labels in a separate plot
            ax=fig.add_subplot(n_params,n_params,label_loc)
            ax.set_xlim([-1,1])
            ax.set_ylim([-1,1])
            
            for i in range(len(fishermat)) :
                ax.plot([-1,1],[-3,-3],color=lc[i],linestyle=ls[i],
                        linewidth=lw[i],label=labels[i])
            ax.legend(loc='upper left',frameon=False,fontsize=FS)
            ax.axis('off')

    if fname!="none" :
        plt.savefig(fname,bbox_inches='tight')

    plt.show()
