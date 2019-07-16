#! /usr/bin/env python
##plotting exercise
##created by Nate
##06/24/2019
import numpy as np
import matplotlib.pyplot as plt

x =np.linspace(0,10,1000)
v1 = x/(1+x)
v1_5 = x/((1+x**(3/2))**(1/(3/2)))
v2 = x/((1+x**2)**(1/2))
v4 = x/((1+x**4)**(1/4))
v10 = x/((1+x**10)**(1/10))
v25 = x/((1+x**25)**(1/25))
v100 = x/((1+x**100)**(1/100))
ve = 1-np.exp(-x)

plt.plot(x, v1, '-r', label='m=1')
#plt.plot(x, v1_5, '-m', label='m=1.5')
plt.plot(x, v2, '-k', label='m=2')
plt.plot(x, v4, '-b', label='m=4')
#plt.plot(x, v10, '-g', label='m=10')
#plt.plot(x, v25, '-c', label='m=25')
plt.plot(x, v100, '-g', label='m=100')
plt.plot(x, ve, '--m', label='exp')

plt.title('Model Rotation Curves')
plt.xlabel('$R/R_O$', color='#1C2833')
plt.ylabel('$V/V_O$', color='#1C2833')
plt.legend(loc='lower right')
plt.grid()
plt.savefig('model_rotcur.png')
plt.show()

