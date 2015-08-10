import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import argparse
import os
import model_builder as mdb

def load_file(cwd,subdir,iteration,temp,filename):
    print 'loading epsilons from iteration %s'%(iteration)
    file_location = '%s/iteration_%s/%s_0'%(subdir,iteration,temp)
    os.chdir(file_location)
    model_params = np.loadtxt(filename)
    os.chdir(cwd)
    return model_params

def filter_params(param_x):
    print 'filtering the parameters'
    delta = []
    try:
        for z in xrange(args.iterations+1):
            delta.append(abs(1-param_x[z]))
            change = max(delta)
        if change < args.min_change:
            param_x = []
    except IndexError:
        pass
    filtered = 'FILTERED'
    return param_x,filtered

def filter_above_limit():
    print 'filtering the parameters that do not reach 0.01'
    x = min(param_x)
    if x > 0.01:
        param_x = []
    else:
        pass
    filtered = 'LIMIT'
    return param_x,filtered

def get_args():
    parser = argparse.ArgumentParser(description = "enter the path, file name; will then create graph of epsilon's progression throughout each iteration")
    parser.add_argument('--subdir')
    parser.add_argument('--file',default='model_params',help='default is model_prams')
    parser.add_argument('--iterations',type=int,help='number of iterations')
    parser.add_argument('--temp',type=int)

    action_group = parser.add_mutually_exclusive_group()
    action_group.add_argument('--trend','-t',action='store_true',default=False)
    action_group.add_argument('--limit','-l',action='store_true',default=False)

    plot_group = parser.add_mutually_exclusive_group()
    plot_group.add_argument('--all','-a',action='store_true',default=False)
    plot_group.add_argument('--one','-o',type=int)
    plot_group.add_argument('--range','-r',nargs=2,type=int)

    parser.add_argument('--filter_switch',default=False,action='store_true')
    parser.add_argument('--min_change',type=float,default=None)
        
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    model,fitopts = mdb.inputs.load_model(args.subdir,dry_run=False)
    cwd = os.getcwd()

    if args.all:
        os.chdir('%s/iteration_0/%s_0'%(args.subdir,args.temp))
        params = np.loadtxt('model_params')
        num_params = [0, params.size]
        os.chdir(cwd)
    elif args.one:
        num_params = [args.one, args.one+1]
    elif args.range:
        num_params = [args.range]
    
    for x in xrange(num_params[0],num_params[1]):
        param_x = [1]
        for i in xrange(args.iterations):
            model_params = load_file(cwd,args.subdir,i,args.temp,args.file)
            param_x.append(model_params[x])
            if args.filter_switch:
                param_x,filtered = filter_params(param_x)
            elif args.limit:
                filter_above_limit()
            else:
                filtered = 'UNFILTERED'
        plt.title('%s_%s_Contacts_%s'%(filtered,args.min_change,num_params))
        plt.plot(param_x)
    plt.savefig('%s_%s_Contacts_%s.png'%(filtered,args.min_change,num_params))
            
    
exit()
