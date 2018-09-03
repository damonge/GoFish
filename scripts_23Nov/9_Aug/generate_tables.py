import functions as f
import numpy as np
 
keywords = ["z0", "norm", "zlo", "zhi", "n_bins", "sigma_z_scale", "galaxy_bias_scale", "nz_lines", "marg_bz", "marg_sz"];
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

# Generate bins_c.txt (redshift bins)
fp = open("bins_c.txt", "w+");
fp.write("# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax")
(lo, hi) = f.generate_bins(n_bins, zlo, zhi, z0, norm);
for i in range(0, n_bins):
    fp.write("\n%f %f %f %d %d %d" % (lo[i], hi[i], f.sigma_z(f.median(lo[i], hi[i], z0), sigma_z_scale), 0, 0, f.lmax(f.median(lo[i], hi[i], z0)) ) );
fp.close();

# Generate bins_s.txt (redshift bins)
fp = open("bins_s.txt", "w+");
fp.write("# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax")
(lo, hi) = f.generate_bins(n_bins, zlo, zhi, z0, norm);
for i in range(0, n_bins):
    fp.write("\n%f %f %f %d %d %d" % (lo[i], hi[i], f.sigma_z(f.median(lo[i], hi[i], z0), sigma_z_scale), 0, 0, 5000 ) );
fp.close();

# Generate nz_c.txt (galaxy distribution)
fp = open("nz_c.txt", "w+");
stepsize = 2*(zhi-zlo)/(nz_lines-1);
for i in range (-int(np.floor(zlo/stepsize)), nz_lines):
    fp.write("%f %f" % (zlo+i*stepsize, f.nz(zlo+i*stepsize, z0, norm)));
    if (i != (nz_lines-1)):
        fp.write("\n");
fp.close();

# Generate bz_c.txt (galaxy bias)
fp = open("bz_c.txt", "w+");
# For some reason z=0 is required, but we don't want to marginalize over this
fp.write("%f %f %d -1\n" % (0, f.galaxy_bias(0, galaxy_bias_scale), 0 ) );
for i in range(0, n_bins):
    fp.write("%f %f %d -1" % (f.median(lo[i], hi[i], z0), f.galaxy_bias(f.median(lo[i], hi[i], z0), galaxy_bias_scale), marg_bz ) );
    if (i != (n_bins - 1)):
        fp.write("\n");
fp.close();

# Generate sz_c.txt (magnification bias)
fp = open("sz_c.txt", "w+");
# For some reason z=0 is required, but we don't want to marginalize over this
fp.write("%f %f %d -1\n" % (0, f.magnification_bias(), 0 ) );
for i in range (0, n_bins):
    fp.write("%f %f %d -1" % (f.median(lo[i], hi[i], z0), f.magnification_bias(), marg_sz) );
    if (i != (n_bins - 1)):
        fp.write("\n");
fp.close();
