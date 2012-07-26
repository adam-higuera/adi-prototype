#!/usr/bin/python

import scipy
import pylab

for i in range(0,100):
    fields=scipy.genfromtxt("dump" + str(i) + ".txt", delimiter=",")
    pylab.hot()
    pylab.pcolor(fields, vmin=-1.0, vmax=1.0);
    pylab.savefig("dump" + str(i) + ".png")
