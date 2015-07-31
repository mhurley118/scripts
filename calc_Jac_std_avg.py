import numpy as np
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('temp')
parser.add_argument('--brs',default=1000,type=int)
args = parser.parse_args()

print 'Loading file of all Jacobians'
Jac = np.loadtxt('Jac_T'+args.temp+'_0'+str(args.brs)+'_resamples.txt')
Jac.shape = (args.brs+1,-1,977)
x,y,z = Jac.shape

print 'Calculating standard deviation'
stddevs = []
for i in xrange(y):
    for j in xrange(z):
        val = std(Jac[:,i,j])
        stddevs.append(val)

print 'Calculatin averages'
averages = []
for i in xrange(y):
    for j in xrange(z):
        val = mean(Jac[:,i,j])
        averages.append(val)

print 'Saving standard deviations to file'
sfile = open('std_T160_0'+str(args.brs)+'_resamples.txt','w+')
for i in xrange(len(stddevs)):
    sfile.write('%s '%(stddevs[i]))

print 'Saving averages to file'
afile = open('avgJac_T160_0'+str(args.brs)+'_resamples.txt','w+')
for i in xrange(len(averages)):
    afile.write('%s '%(averages[i]))

print 'plotting standard deviations'
try:
    os.mkdir('plots')
except:
    pass
os.chdir('plots')
stddevs = np.asarray(stddevs)
stddevs.shape = (-1, 977)
print stddevs.shape
fig = plt.figure()
ax = fig.add_subplot(111, title='Standard Deviations')
im = ax.pcolor(stddevs, cmap = 'Blues')
fig.colorbar(im)
fig.savefig('stddev')

