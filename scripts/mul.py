#! /usr/bin/env python
#
#  example access to a table, meant for subsetting using boolean expressions on columns
#

import numpy as np


# method 1: loadtxt().T ============================
# i.e. how the table looks in a text file

data = np.loadtxt('mul.tab').T
# array([[ 1. ,  1. ,  1. ,  1. ,  2. ,  2. ,  2. ,  2. ],
#        [ 1. ,  2. ,  3. ,  4. ,  1.1,  2.2,  3.3,  4.4],
#        [ 2. ,  3.1,  2.9,  4.2,  4.4,  4.9,  5.4,  6.1]])

m1 = np.where(data[0] == 1)
m2 = np.where(data[0] == 2)
# m1 -> (array([0, 1, 2, 3]),)






# method 2: loadtxt() ===============================

# 1) get data
data = np.loadtxt('mul.tab')
# array([[ 1. ,  1. ,  2. ],
#        [ 1. ,  2. ,  3.1],
#        [ 1. ,  3. ,  2.9],
#        [ 1. ,  4. ,  4.2],
#        [ 2. ,  1.1,  4.4],
#        [ 2. ,  2.2,  4.9],
#        [ 2. ,  3.3,  5.4],
#        [ 2. ,  4.4,  6.1]])

# 2) get index
m1 = np.where(data[:,0] == 1)   
m2 = np.where(data[:,0] == 2)
# m1 -> (array([0, 1, 2, 3]),)

# 3) subset data 
data1 = data[m1]
data2 = data[m2]
# array([[ 1. ,  1. ,  2. ],
#        [ 1. ,  2. ,  3.1],
#        [ 1. ,  3. ,  2.9],
#        [ 1. ,  4. ,  4.2]])

# 4) create vectors we want for plotting
x1 = data1[1]
y1 = data1[2]
s1 = y1/s1

x2 = data2[1]
y2 = data2[2]

# in python 2) 3) and 4) this could be combined into one method:
# x1 = data[np.where(data[:,0] == 1)][1]
# y1 = data[np.where(data[:,0] == 1)][2]




# method3 : np.extract =====================================

# 1) get data
data = np.loadtxt('mul.tab').T
m = data[0]
x = data[1]
y = data[2]

or better:

# 1) get data
(m,x,y) = np.loadtxt('mul.tab').T

# 2) get vectors for plotting
x1 = np.extract(m==1, x)
y1 = np.extract(m==1, y)
s1 = np.extract(m==1, y/x)

x2 = np.extract(m==2, x)
y2 = np.extract(m==2, y)


x3 = np.extract((m==1) | (m==2), x)
y3 = np.extract((m==1) | (m==2), y)
