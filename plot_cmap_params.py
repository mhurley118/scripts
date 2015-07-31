""" Plot model parameters for a given iteration """

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import os

def plot_color(pairs_U, pairs_L, params_U, params_L, iter):
    plt.figure()
    x_U = pairs_U[:,0]
    y_U = pairs_U[:,1]
    z_U = params_U
    
    x_L = pairs_L[:,0]
    y_L = pairs_L[:,1]
    z_L = params_L
    
    if np.min(z_U) < np.min(z_L):
        zmin = np.min(z_U)
    else:
        zmin = np.min(z_L)
    
    z_mean = np.mean(z_L)
    zmax = z_mean + (z_mean - zmin)
    
    # Set constant zmin, zmax
    zmin = 0.2
    zmax = 1.8
    
    cp = plt.scatter(x_U, y_U, s=8, c=z_U, marker='o', linewidth=0., vmin=zmin, vmax=zmax)
    cb = plt.colorbar(cp)
#    cb.set_label("Contact Epsilons")
    cp = plt.scatter(y_L, x_L, s=8, c=z_L, marker='o', linewidth=0., vmin=zmin, vmax=zmax)
    maxval = np.max(x_U)
    plt.axis([0, maxval, 0, maxval])
    plt.xlabel("Residue i")
    plt.ylabel("Residue j")
    plt.title("Iteration %d Epsilons"%iter)
    plt.savefig("epsilons_iter_%d.png"%iter)

    plt.close()
    
def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('--upairs', default='pairwise_params', type=str, help='Pairs to plot in upper triangle')
    parser.add_argument('--lpairs', default='pairwise_params', type=str, help='Pairs to plot in lower triangle')
    parser.add_argument('--uparams', default='model_params', type=str, help='Params to plot in upper triangle')
    parser.add_argument('--lparams', default=None, type=str, help='Params to plot in lower triangle')
    parser.add_argument('--subdir', type=str, help='Directory containing files')
    parser.add_argument('--iter', type=int, help='Iteration')

    args = parser.parse_args()

    return args
    
if __name__=="__main__":

    args = get_args()
    
    iter = args.iter
    
    location = "%s/iteration_%d"%(args.subdir, iter)
    os.chdir(location)
    
    pairs_U = np.loadtxt("%s"%args.upairs, dtype=int, usecols=(0,1))
    pairs_L = np.loadtxt("%s"%args.lpairs, dtype=int, usecols=(0,1))
    params_U = np.loadtxt("%s"%args.uparams, dtype=float, usecols=(0,))
    
    if args.lparams == None:
        params_L = np.ones(np.shape(pairs_L)[0])
    else:
        params_L = np.loadtxt("%s"%args.lparams, dtype=float, usecols=(0,))

    plot_color(pairs_U, pairs_L, params_U, params_L, iter)
