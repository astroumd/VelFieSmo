#! /usr/bin/env python
#
##plotting exercise
##created by Nate
##06/24/2019
import numpy as np
import matplotlib.pyplot as plt

def vrot(x, m):
    return x / pow(1+x**m, 1/m)

x =np.linspace(0,10,1000)

ve = 1-np.exp(-x)

v7 = vrot(x,0.75)
v1 = vrot(x,1)
v2 = vrot(x,2)
v4 = vrot(x,4)
v100 = vrot(x,100)

#plt.plot(x, v7, '-r', label='m=0.75')
plt.plot(x, v1, '-r', label='m=1')
#plt.plot(x, v1_5, '-m', label='m=1.5')
plt.plot(x, v2, '-k', label='m=2')
plt.plot(x, v4, '-b', label='m=4')
#plt.plot(x, v10, '-g', label='m=10')
#plt.plot(x, v25, '-c', label='m=25')
plt.plot(x, v100, '-g', label='m=100')
plt.plot(x, ve, '--m', label='exp')

plt.title('Model Rotation Curves: $V/V_0 = r / (1+r^m)^{1/m},   r=R/R_0$')
plt.xlabel('$R/R_O$')
plt.ylabel('$V/V_O$') 
plt.legend(loc='lower right')
plt.grid()
plt.savefig('model_rot_curve.pdf')
plt.show()

