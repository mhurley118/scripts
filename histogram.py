import argparse
from numpy import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()
parser.add_argument('file')
parser.add_argument('--path', default=None, help='if file is not in same directory')
parser.add_argument('--bins', default=100, type=int)
parser.add_argument('title')
args = parser.parse_args()

print os.getcwd()
x_val = []
print 'determining x values'
for i in xrange(args.bins):
    x_val.append(i*0.1)

print 'loading data'
data = np.loadtxt(args.file)
data = data.ravel()
data = data.tolist()
#print data

print 'binning data'
count = []
for q in xrange(args.bins):
    list_q = []
    for i in xrange(len(data)):
        a = 0.1*q
        b = 0.1*(q+1)
        if data[i] < b and data[i] >= a:
            list_q.append(data[i])
#    print 'length='+str(len(list_q))
    count.append(len(list_q))

os.chdir('/home/mjh9/scratch/output_files')
fig = plt.figure()
ax = fig.add_subplot(111)
for i in xrange(args.bins):
    ax.bar(x_val[i], count[i], width=0.1)
fig.savefig(args.title+'.png')

exit()
