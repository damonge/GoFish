import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
from scipy.interpolate import interp1d as ip
from astropy.cosmology import Planck15 as cosmo

def integral_nz_unnormalized(zmin, zmax, z0):
    hi = (2.*z0/3.)*np.exp(-(zmax/z0)**1.5)*(-(zmax/z0)**1.5-1.);
    lo = (2.*z0/3.)*np.exp(-(zmin/z0)**1.5)*(-(zmin/z0)**1.5-1.);
    return (hi - lo);

def normalization(norm, z0):
    return norm/integral_nz_unnormalized(0., 100., z0);

def nz(z, z0, norm):
    return normalization(norm, z0) * (z/z0)**2 * np.exp(-(z/z0)**1.5);

def integral_nz(zmin, zmax, z0, norm):
    return normalization(norm, z0)*integral_nz_unnormalized(zmin, zmax, z0);

def solve_for_hi(lo, norm, n_bins, zlo, zhi, z0):
    area = integral_nz(zlo, zhi, z0, norm)/n_bins;
    def reduced_integral_nz(hi):
        return abs(integral_nz(lo, hi, z0, norm) - area);
    hi = opt.minimize_scalar(reduced_integral_nz, method = 'bounded', bounds=(lo, zhi)).x;
    return hi;

def median(zlo, zhi, z0):
    def area_difference(z):
        return abs(integral_nz_unnormalized(zlo, z, z0) - integral_nz_unnormalized(z, zhi, z0));
    z = opt.minimize_scalar(area_difference, method = 'bounded', bounds=(zlo, zhi)).x;
    return z;

def generate_bins(n_bins, zlo, zhi, z0, norm):
    lo = np.zeros(n_bins);
    hi = np.zeros(n_bins);
    lo[0] = zlo;
    hi[-1] = zhi;
    for i in range(0, n_bins-1):
        hi[i] = solve_for_hi(lo[i], norm, n_bins, zlo, zhi, z0);
        lo[i+1] = hi[i];
    return (lo, hi);

def sigma_z(z, scale):
    return scale*(1. + z);

#def chi(z):
#    return cosmo.comoving_distance(z).value / cosmo.h;

# Import the values (z, lmax) from Joachimi, Bridle 2010
fp = open("lmax_Joachimi.txt", "r");
z_joachimi = [];
l_max_joachimi = [];
import re
for line in fp:
    if re.search('[a-zA-Z]', line):
        continue;
    else:
        first_space = line.find(' ');
        z_joachimi.append(float(line[:first_space]));
        line = line[first_space:];
        while line.find(' ')==0:
            line = line[1:];
        l_max_joachimi.append(float(line));

def lmax(z):
    # Interpolate from Joachimi's values
    interpolation = ip(z_joachimi, l_max_joachimi);
    return(interpolation(z));
    #k_max = 0.132 * z * cosmo.h; # Mpc^-1, Joachimi & Bridle 2010
    #return chi(z)*k_max;

def galaxy_bias(z, scale):
    return scale;

def magnification_bias():
    alpha = 0.7; # fiducial value from Duncan et al 2013
    s = - (2./5.) * alpha; # not sure about this minus sign -- LFT May17
    return s;
