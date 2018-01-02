from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

x = np.linspace(0, 10, num=11, endpoint=True)
y = np.cos(-x**2/9.0)
f = interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')
xnew = np.linspace(0, 10, num=41, endpoint=True)
result = integrate.quad(lambda x:f(x), 0, 4.5)
print result

plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
plt.legend(['data', 'linear', 'cubic'], loc='best')
plt.show()