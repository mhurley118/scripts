from numpy import *
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import mdtraj as md
import argparse
from time import strftime
import os
import model_builder

def bootstrap(temp,fitopts):
    print 'configuring directories'
    cwd = args.cwd
    subdir = args.subdir
    iteration = fitopts['iteration']
    traj_location = '%s/iteration_%s/%s_0'%(subdir,iteration,temp)
    os.chdir(traj_location)
    
    print 'calculating distances'
    print os.getcwd()
    FRET = md.load('traj.xtc', top = 'Native.pdb')
    oFRET = md.compute_distances(FRET, np.asarray(args.pairs),periodic=False)
    oFRET = oFRET.ravel()

    print 'binning data'
    count_one =[]
    for q in xrange(args.bins):
        list_q = []
        for i in xrange(oFRET.size):
            a = 0.1*q
            b = 0.1*(q+1)
            if oFRET[i] < b and oFRET[i] >= a:
                list_q.append(oFRET[i])
        count_one.append(len(list_q))

    print 'bootstrapping'
    for i in xrange(args.resamplings):
        print 'RESAMPLE: '
        print i
        #this loop controls how many iterations are run--5 bootstraps, 10, or 10,000 bootstraps
        count_two = []
        rFRET = npr.choice(oFRET, oFRET.size)
        for q in xrange(args.bins):
            #the following portion is the histogram code; q essentially controls the number of bins
            list_q = []
            for i in xrange(len(rFRET)):
                #this runs all 100,001 numbers through a filter 0.1 nm in width
                a = 0.1*q
                b = 0.1*(q+1)
                if rFRET[i] < b and rFRET[i] >= a:
                    list_q.append(rFRET[i])
            count_two.append(len(list_q))
        for j in xrange(len(count_two)):
            count_one.append(count_two[j])
    print 'setting up data into array'
    count_one = np.asarray(count_one)
    print 'shaping array'
    count_one.shape = (-1,args.bins)
    os.chdir(cwd)
    return count_one

def plot_and_save(count,temp,model_type):
    cwd = args.cwd
    #the average number of counts in each bin
    print 'finding the avg count in each bin'
    avg_counts = []
    for i in xrange(args.bins):
        avg_counts.append(mean(count[:,i]))
    print 'calculating std'
    stddevs = []
    for i in xrange(args.bins):
        #find the standard deviations
        stddevs.append(std(count[:,i]))
    
    date = strftime('_%b_%d_')

    print 'creating figure'
    x_val = []
    print 'determining x values'
    for i in xrange(args.bins):
        x_val.append(i*0.1)
    
    os.chdir(args.savelocation)
    try:
        os.mkdir('mdtraj_dist'+date+'T')
    except:
        pass
    os.chdir('mdtraj_dist'+date+'T')
    fig = plt.figure()
    ax = fig.add_subplot(111,title='Histogram of Avg Distances between '+str(args.pairs)+' at T'+str(temp)+','+str(args.resamplings)+' Resamplings',ylabel='Count',xlabel='Nanometers',ylim=[0,max(avg_counts)+100])
    print max(avg_counts)
    for i in xrange(args.bins):
        ax.bar(x_val[i], avg_counts[i], width = 0.1, yerr = stddevs[i], color = 'limegreen', ecolor = 'k')
    fig.savefig('%s,T%s_%s_avg_dist.png'%(model_type,temp,args.resamplings))

    print 'creating files'
    np.savetxt('%s,T%s_%s_avg_counts'%(model_type,temp,args.resamplings),avg_counts)
    np.savetxt('%s,T%s_%s_stddevs'%(model_type,temp,args.resamplings), stddevs)
    os.chdir(cwd)

def sanitize_args(args):
    original_directory = os.getcwd()
    os.chdir(args.cwd)
    ##If temps was not specified, will attempt to open Temparray.txt, Otherwise prints ERROR
    if args.temps==None and os.path.isfile("Temparray.txt"):
        args.temps = np.loadtxt("Temparray.txt",dtype=int)
    else:
        print "ERROR: No Temperature Directories Specified"
    
    ##set the pairs for fitting in the array format for the calculation
    pairs = np.array([[args.pairs[0], args.pairs[1]]])
    print "Number of pairs is: %d" % len(args.pairs)
    if len(args.pairs)>2:
        for i in np.arange(3, len(args.pairs), 2):
            pairs = np.append(pairs, np.array([[args.pairs[i-1], args.pairs[i]]]), axis=0)
    args.pairs = pairs
        
    os.chdir(original_directory)
    return args

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subdir',default='1PB7',help='Subdirectory, also name of .ini file')
    parser.add_argument('--temps',default=None,type=int,nargs='+')
    parser.add_argument('--cwd',default=os.getcwd(),help='directory for fitting')
    parser.add_argument('--pairs',nargs='+',type=int,default=[219, 371],help='FRET pairs')
    parser.add_argument('--bins',type=int,default=100,help='Number of bins in histogram')
    parser.add_argument('--resamplings','--resamples',type=int,default=1000,help='Number of resamplings')
    parser.add_argument('--savelocation',default=os.getcwd())
    parser.add_argument('--model_type',default='CACB')

    args = parser.parse_args()
    args = sanitize_args(args)
    return args

if __name__ == '__main__':
    args = get_args()
    model,fitopts = model_builder.inputs.load_model(args.subdir,dry_run=False)
    
    for i in xrange(len(args.temps)):
        count = bootstrap(args.temps[i],fitopts)
        plot_and_save(count,args.temps[i],args.model_type)

exit()
