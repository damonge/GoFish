import functions as f
import numpy as np
 
keywords = ["z0", "norm", "zlo", "zhi", "n_bins", "sigma_z_scale", "galaxy_bias_scale", "nz_lines", "marg_bz", "marg_sz", "nbins_shear", "nbins_clustering"];
parameters = np.zeros(len(keywords));
for i in range(0, len(keywords)):
    fp = open("input.txt", "r");
    string = keywords[i];
    for line in fp:
        if string in line:
            snippet = line[len(string)+1:];
            end = snippet.find(" ");
            if (end != -1):
                parameters[i] = float(snippet[0:end]);
            else:
                parameters[i] = float(snippet);

z0 = parameters[0];
norm = parameters[1]; #arcmin^-2
zlo = parameters[2]; #boundaries of the z under consideration
zhi = parameters[3];
n_bins = int(parameters[4]);
sigma_z_scale = parameters[5];
galaxy_bias_scale = parameters[6];
nz_lines = int(parameters[7]);
marg_bz = int(parameters[8]);
marg_sz = int(parameters[9]);
n_bins_shear = int(parameters[10]);
n_bins_clustering = int(parameters[11]);

# Generate bins_c.txt (redshift bins)
#(lo, hi) = f.generate_bins(n_bins, zlo, zhi, z0, norm);
#for i in range(0, n_bins):
#    fp.write("\n%f %f %f %d %d %d" % (lo[i], hi[i], f.sigma_z(f.median(lo[i], hi[i], z0), sigma_z_scale), 0, 0, f.lmax(f.median(lo[i], hi[i], z0)) ) );
fp_bin_boundaries = np.loadtxt('bin_boundaries.txt')
shear_lines = fp_bin_boundaries[:n_bins_shear,:]
for ii in xrange(n_bins_shear):
    fp = open("bins_s_"+str(ii)+".txt", "w+");
    fp.write("# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax")
    #fp.write("\n%f %f %f %d %d %d" % (shear_lines[ii,1], shear_lines[ii,2], sigma_z_scale, 0, 0, 5000))
    fp.write("\n%f %f %f %d %d %d" % (0., 3., sigma_z_scale, 0, 0, 5000))
    fp.close();
clustering_lines = fp_bin_boundaries[n_bins_shear:n_bins_shear+n_bins_clustering,:]
for ii in xrange(n_bins_clustering):
    fp = open("bins_c_"+str(ii)+".txt", "w+");
    fp.write("# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax")
    #fp.write("\n%f %f %f %d %d %d" % (clustering_lines[ii,1], clustering_lines[ii,2], sigma_z_scale, 0, 0, f.lmax(clustering_lines[ii,3])))
    # TODO : this is just a temporary fix to neutralize the potential bug in CLASS
    fp.write("\n%f %f %f %d %d %d" % (0., 3., sigma_z_scale, 0, 0, f.lmax(clustering_lines[ii,3])))
    fp.close();

# Generate bins_s.txt (redshift bins)
#fp = open("bins_s.txt", "w+");
#fp.write("# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax")
#(lo, hi) = f.generate_bins(n_bins, zlo, zhi, z0, norm);
#for i in range(0, n_bins):
#    fp.write("\n%f %f %f %d %d %d" % (lo[i], hi[i], f.sigma_z(f.median(lo[i], hi[i], z0), sigma_z_scale), 0, 0, 5000 ) );
#fp.close();

# Generate nz_c.txt (galaxy distribution)
#fp = open("nz_c.txt", "w+");
#stepsize = 2*(zhi-zlo)/(nz_lines-1);
#for i in range (-int(np.floor(zlo/stepsize)), nz_lines):
#    fp.write("%f %f" % (zlo+i*stepsize, f.nz(zlo+i*stepsize, z0, norm)));
#    if (i != (nz_lines-1)):
#        fp.write("\n");
#fp.close();


# generate source number densities
DES_nz_file = np.load('DES_nz.npz')
l = DES_nz_file['source']
z_centres = l['Z_MID']
for ii in xrange(n_bins_shear):
   # TODO : temporary replacement 
   if ii == 0 :
      fp = open('nz_shear_0_modified.txt', 'w+')
      fp.write("%f %f\n" % (0.,0.))
      for jj in xrange(len(z_centres)):
         fp.write("%f %f\n" % (z_centres[jj],5.*l['BIN'+str(ii+1)][jj]))
      fp.close()
   if ii == 1 :
      fp = open('nz_shear_1_modified.txt', 'w+')
      fp.write("%f %f\n" % (0.,0.))
      for jj in xrange(len(z_centres)):
         fp.write("%f %f\n" % (z_centres[jj],5.*l['BIN'+str(ii+1)][jj]))
      fp.close()

   fp = open('nz_shear_' + str(ii) + '.txt', 'w+')
   fp.write("%f %f\n" % (0.,0.))
   for jj in xrange(len(z_centres)):
      fp.write("%f %f\n" % (z_centres[jj],l['BIN'+str(ii+1)][jj]))
   fp.close()
# generate lens number densities
DES_nz_file = np.load('DES_nz.npz')
l = DES_nz_file['lens']
z_centres = l['Z_MID']
for ii in xrange(n_bins_clustering):
   # TODO : temporary replacement 
   if ii == 0 :
      fp = open('nz_clustering_0_modified.txt', 'w+')
      fp.write("%f %f\n" % (0.,0.))
      for jj in xrange(len(z_centres)):
         fp.write("%f %f\n" % (z_centres[jj],5.*l['BIN'+str(ii+1)][jj]))
      fp.close()
   if ii == 1 :
      fp = open('nz_clustering_1_modified.txt', 'w+')
      fp.write("%f %f\n" % (0.,0.))
      for jj in xrange(len(z_centres)):
         fp.write("%f %f\n" % (z_centres[jj],5.*l['BIN'+str(ii+1)][jj]))
      fp.close()

   fp = open('nz_clustering_' + str(ii) + '.txt', 'w+')
   fp.write("%f %f\n" % (0.,0.))
   for jj in xrange(len(z_centres)):
      fp.write("%f %f\n" % (z_centres[jj],l['BIN'+str(ii+1)][jj]))
   fp.close()

# Generate bz_c.txt (galaxy bias)
fp = open("bz_c.txt", "w+");
# For some reason z=0 is required, but we don't want to marginalize over this
fp.write("%f %f %d -1\n" % (0, f.galaxy_bias(0, galaxy_bias_scale), 0 ) );
for ii in range(0, n_bins_clustering):
    fp.write("%f %f %d -1" % (clustering_lines[ii,3], f.galaxy_bias(clustering_lines[ii,3], galaxy_bias_scale), marg_bz ) );
    if (i != (n_bins - 1)):
        fp.write("\n");
fp.close();

# Generate sz_c.txt (magnification bias)
fp = open("sz_c.txt", "w+");
fp_alpha_values = np.loadtxt('alpha_values.txt')
# For some reason z=0 is required, but we don't want to marginalize over this
fp.write("%f %f %d -1\n" % (0, 0.4*fp_alpha_values[0,1], 0 ) );
for ii in range (0, n_bins_clustering):
    fp.write("%f %f %d -1" % (clustering_lines[ii,3], 0.4*fp_alpha_values[ii,1], marg_sz) );
    if (ii != (n_bins - 1)):
        fp.write("\n");
fp.close();
