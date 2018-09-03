import numpy as np
import matplotlib.pyplot as plt

dndn_cases = ['G']

edn_cases = ['d']

nrows_cases = ['0', '1', '2', '3', '4', '5']


# om : 0.316
# s8 : 0.831
# ob : 0.049
# hh : 0.690
# ns : 0.966
# ok : 0.000
# w0 : -1.000

params = [
    '$\Omega_M$',
    '$\sigma_8$',
    '$\Omega_B$',
    '$h$',
    '$n_s$',
    '$\Omega_K$',
    '$w_0$'
    ]

cases = [
    'd0',
    'd1',
    'd2',
    'd3',
    'd4',
    'd5'
        ]

biases_A = np.zeros((len(params), len(cases)))
errors_A = np.zeros((len(params), len(cases)))
ratios_A = np.zeros((len(params), len(cases)))

indices = np.linspace(1., float(len(cases)), num = len(cases))

nn = 0
# loop over the given cutting schemes
for ii in range(0, len(dndn_cases)):
    for jj in range(0, len(edn_cases)):
        for kk in range(0, len(nrows_cases)):
            folder = dndn_cases[ii] + '_' + edn_cases[jj] + '_' + nrows_cases[kk] + '_no/'
            fisher_object = np.load(folder + 'fisher_raw.npz')
            bias_vec = fisher_object['fisher_bias']
            fisher_matrix = fisher_object['fisher_tot']
            inverse = np.linalg.inv(fisher_matrix)
            biases = abs(np.dot(inverse, bias_vec))[:len(params)]
            print 1e2*biases
            errors = np.zeros(len(biases))[:len(params)]
            ratios = np.zeros(len(biases))[:len(params)]
            for mm in range(0, len(biases)):
                errors[mm] = np.sqrt(inverse[mm,mm])
                ratios[mm] = biases[mm]/errors[mm]
            biases_A[:,nn] = biases
            errors_A[:,nn] = errors
            ratios_A[:,nn] = ratios
            nn += 1

TITLE_SIZE = 16
LABEL_SIZE = 16
plt.rc('axes', titlesize=TITLE_SIZE)
plt.rc('axes', labelsize=TITLE_SIZE)

plt.figure(figsize=(30,15))
for index in range(0, len(params)):
    plt.subplot(3, len(params), 1 + index)
    plt.plot(indices, 1e2*biases_A[index,:])
    plt.xticks(indices, cases)
    if index == 0:
        plt.ylabel('abs(bias) (*1e2)')
    plt.title(params[index])
    plt.subplot(3, len(params), 1+len(params) + index)
    plt.plot(indices, 1e2*errors_A[index,:])
    plt.xticks(indices, cases)
    if index == 0:
        plt.ylabel('error (*1e2)')
    plt.subplot(3, len(params), 1+2*len(params) + index)
    plt.plot(indices, ratios_A[index,:], label = 'G')
    plt.xticks(indices, cases)
    if index == 0:
        plt.ylabel('ratio bias/error')
    if index == len(params)-1:
        plt.legend(loc = 'lower right')


#plt.show()


plt.savefig('figure1.pdf', bbox_inches='tight')

