#! /usr/bin/env python
#
#   effect  of beam smearing on inclination

import numpy as np
import matplotlib.pyplot as plt


# data are taken from rotcurshape on variations on tghe following
# rm -rf model0; ../scripts/mkgalcube m=2 beam=1 inc=75 re=1 vel=mom

i1  = np.array([60,   65,   70,   75,   80])
i1f = np.array([58.3, 63.0, 67.5, 71.7, 75.1])
di1 = i1-i1f

# re=1000
i2  = np.array([75])
i2f = np.array([74.1])
di2 = i2-i2f

# re=1
i3  = np.array([80,   75,   70,   65,   60])
i3f = np.array([72.1, 69.5, 65.8, 61.6, 57.2])
di3 = i3-i3f

# re=0.5
i4  = np.array([75])
i4f = np.array([66.3])
di4 = i4-i4f


plt.plot(i1,       di1, label='re=2')
plt.scatter(i2,    di2, label='re=1000')
plt.scatter(i3,    di3, label='re=1')
plt.scatter(i4,    di4, label='re=0.5')

plt.title('Kinematic Inclination; m=2')
plt.xlabel('INC')
plt.ylabel('INC-FIT')
plt.legend(loc='upper left')
plt.grid()
plt.savefig('smooth_inc.png')
plt.show()

