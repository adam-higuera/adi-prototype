#!/usr/local/bin/python

import scipy
import pylab

fields=scipy.genfromtxt("dump3.txt", delimiter=",")
pylab.hot()
pylab.pcolor(fields, vmin=-1.0, vmax=1.0);
pylab.show()
