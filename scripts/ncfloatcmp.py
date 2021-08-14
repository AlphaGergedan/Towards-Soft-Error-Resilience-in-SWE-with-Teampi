##--------------------------------------------------##
##  nc file comparison script for debugging         ##
##  Compares the arrays of two given SWE solutions  ##
##  usage: python ncfloatcmp.py file1.nc file2.nc   ##
##  author: Atamert Rahma | rahma@in.tum.de         ##
##--------------------------------------------------##

import numpy as np
import netCDF4 as nc
import sys
import math

fileName1= ""
fileName1 = sys.argv[1]

fileName2= ""
fileName2 = sys.argv[2]

f1 = nc.Dataset(fileName1, 'r')
f2 = nc.Dataset(fileName2, 'r')

t1 = f1.variables['time']
x1 = f1.variables['x']
y1 = f1.variables['y']
b1 = f1.variables['b']
h1 = f1.variables['h']
hu1 = f1.variables['hu']
hv1 = f1.variables['hv']

t2 = f1.variables['time']
x2 = f1.variables['x']
y2 = f1.variables['y']
b2 = f2.variables['b']
h2 = f2.variables['h']
hu2 = f2.variables['hu']
hv2 = f2.variables['hv']

totalfails = 0
failed_y = []
failed_x = []
failed_timeStep = []
failed_array = []

# assert t sizes
assert(len(t1) == len(t2))
t_len = len(t1)
# assert b,h,hu,hv sizes
assert(len(h1) == t_len)
assert(len(h2) == t_len)
assert(len(hu1) == t_len)
assert(len(hu2) == t_len)
assert(len(hv1) == t_len)
assert(len(hv2) == t_len)

# assert x,y sizes
assert(len(x1) == len(x2))
assert(len(y1) == len(y2))
x_len = len(x1)
y_len = len(y1)
assert(len(b1) == y_len)
assert(len(b1) == y_len)
################################
# assert h,hu,hv arrays equal! #
################################
for i in range(0, t_len):
    assert(len(h1[i]) == y_len)
    assert(len(h2[i]) == y_len)
    assert(len(hu1[i]) == y_len)
    assert(len(hu2[i]) == y_len)
    assert(len(hv1[i]) == y_len)
    assert(len(hv2[i]) == y_len)
    for y in range(0, y_len):
        assert(len(h1[i][y]) == x_len)
        assert(len(h2[i][y]) == x_len)
        assert(len(hu1[i][y]) == x_len)
        assert(len(hu2[i][y]) == x_len)
        assert(len(hv1[i][y]) == x_len)
        assert(len(hv2[i][y]) == x_len)
        for x in range(0, x_len):
            h_equal = h1[i][y][x] == h2[i][y][x]
            hu_equal = hu1[i][y][x] == hu2[i][y][x]
            hv_equal = hv1[i][y][x] == hv2[i][y][x]
            if not h_equal:
                totalfails = totalfails + 1
                print("h[" + str(i) + "][" + str(y) + "][" + str(x) + "] : " + str(h1[i][y][x]) + " != " + str(h2[i][y][x]))
                print("timestep number " + str(i) + " = " + str(t1[i]) + "==" + str(t2[i]))
                failed_x.append(x)
                failed_y.append(y)
                failed_timeStep.append(i)
                failed_array.append('h')
            if not hu_equal:
                totalfails = totalfails + 1
                print("hu[" + str(i) + "][" + str(y) + "][" + str(x) + "] : " + str(hu1[i][y][x]) + " != " + str(hu2[i][y][x]))
                print("timestep number " + str(i) + " = " + str(t1[i]) + "==" + str(t2[i]))
                failed_x.append(x)
                failed_y.append(y)
                failed_timeStep.append(i)
                failed_array.append('hu')
            if not hv_equal:
                totalfails = totalfails + 1
                print("hv[" + str(i) + "][" + str(y) + "][" + str(x) + "] : " + str(hv1[i][y][x]) + " != " + str(hv2[i][y][x]))
                print("timestep number " + str(i) + " = " + str(t1[i]) + "==" + str(t2[i]))
                failed_x.append(x)
                failed_y.append(y)
                failed_timeStep.append(i)
                failed_array.append('hv')

##########################
# assert b arrays equal! #
##########################
for y in range(0, y_len):
    assert(len(b1[y]) == x_len)
    assert(len(b2[y]) == x_len)
    for x in range(0, x_len):
        assert(b1[y][x] == b2[y][x])




print("total failed exact comparisons : " + str(totalfails))
print("=======================================================")
for i in range(0,len(failed_x)):
    print("-->  " + str(failed_array[i]) + "[" + str(failed_y[i]) + "][" + str(failed_x[i]) + "] at simulation time " + str(failed_timeStep[i]))
