from numpy import *
import numpy as np
import numpy.random as npr
import os
import time
import argparse
import mdtraj as md
import scipy.stats as stats
import model_builder as mdb
from project_tools.parameter_fitting.util.util import *

model, fitopts = mdb.inputs.load_model('1PB7', dry_run = True)

parser = argparse.ArgumentParser()
parser.add_argument('trajfile')
parser.add_argument('topfile')
parser.add_argument('tempfile')
args = parser.parse_args()

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621
def_FRET_pairs = [[114,192]]
defspacing = 0.1 ## in nm
cwd = os.getcwd()

brs = 1000

def find_sim_bins(savelocation, FRETr, fit_temp, residues=def_FRET_pairs, spacing=defspacing, weights=None):
    #if not weighted, set weights to ones
    if weights == None:
        weights = np.ones(np.shape(FRETr)[0])
    
    #actually histogram it
    print '*********************'
    print np.shape(FRETr)
    print '*********************'
    hist, edges, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)
    hist = hist/(np.sum(hist)*spacing)
    bincenters = 0.5 * (edges[1:] + edges[:-1])
    
    os.chdir(cwd)
    print "Calculated bins for simulation data at a spacing of %.4f" % spacing
    
    return hist, slices

def calculate_average_Jacobian(model,fitopts, FRET_pairs=def_FRET_pairs, spacing=defspacing ):
    if "t_fit" in fitopts:
        fit_temp = fitopts["t_fit"]
    else:
        raise IOError("Missing the fit_temperature, please specify in .ini file")    
    if "fret_pairs" in fitopts:
        fret_pairs = fitopts["fret_pairs"]
        FRET_pairs = np.array(fret_pairs) - 1
        print "The FRET pairs are:"
        print FRET_pairs
    if "y_shift" in fitopts:
        y_shift = fitopts["y_shift"]   
    else:
        y_shift = 0.0
        fitopts["y_shift"] = 0.0        
    if "spacing" in fitopts:
        spacing = fitopts["spacing"]

    ##Define location of logical files    
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    sub = "%s/%s/iteration_%d/" % (cwd,subdir,iteration)
    traj_location = "%s%s/" %(sub, args.tempfile)
    sim_location = "%s/fitting_%d" % (sub,iteration)
        
    os.chdir(traj_location)
    ## Get trajectory, state indicators, contact energy
    print 'Working on calculating model"s trajectory and contact info'
    traj,rij,qij = get_rij_Vp(model)

    ##no procedure included in this code for adding in a y-shift

    sim_feature, sim_slices = find_sim_bins(sim_location, data[:,0], fit_temp, residues=FRET_pairs, spacing=spacing, weights=None)
    beta = 1.0 / (GAS_CONSTANT_KJ_MOL*float(fit_temp))
    
    os.chdir(cwd)
    print "Computing Jacobian and Simparams for the temperature %d, with spacing %f" % (fit_temp, spacing)
    Jacobian = compute_Jacobian_basic(qij,sim_feature*spacing, sim_slices, beta)
    Jacobian = Jacobian/ spacing        
    sim_feature_err = sim_feature ** 0.5
    Jacobian_err = np.zeros(np.shape(Jacobian))
    
    return Jacobian, Jacobian_err, sim_feature, sim_feature_err

def compute_Jacobian_basic(qij, fr, sim_slices, beta, weights=None):
    nbins = np.shape(fr)[0]
    (N_total_traj, npairs) = np.shape(qij)
    if weights == None:
        N_total_weight = N_total_traj
        weights = np.ones(N_total_traj)
    else:
        if not np.shape(weights)[0] == N_total_traj:
            raise IOError("Not every frame is weighted, aborting! Check to make sure weights is same length as trajectory")
        N_total_weight = np.sum(weights)
    
    Jacobian = np.zeros((nbins, npairs),float)
    for idx, bin_location in enumerate(sim_slices):
        Jacobian[bin_location-1, :] += qij[idx,:]*weights[idx]
    Jacobian /= N_total_weight
    Qavg = np.sum(Jacobian, axis=0)
    
    avg_matrix = np.dot(np.array([fr]).transpose(), np.array([Qavg]))
    print "The shape of these matrices are:"
    print np.shape(avg_matrix)
    print np.shape(Jacobian)
    Jacobian -= avg_matrix
    
    Jacobian *= (-1.0) * beta

    return Jacobian

# load file
if 'fret_pairs' in fitopts:
    fret_pairs = fitopts['fret_pairs']
    FRET_pairs = np.asarray(fret_pairs)-1
    print FRET_pairs
t = md.load(args.trajfile, top = args.topfile)
FRETr = md.compute_distances(t, FRET_pairs, periodic=False)
shape = FRETr.shape

maxval = int((np.amax(FRETr)+1)/defspacing)
minval = int((np.amin(FRETr)-1)/defspacing)
num_bins = int(maxval - minval)
ran_size = (int(minval*defspacing),int(maxval*defspacing))

matrix_to_list = []
originalJac = []
for i in xrange(brs+1):
    print i
    data = FRETr
    if i < 1:
        j, je, sf, sfe = calculate_average_Jacobian(model,fitopts, FRET_pairs=def_FRET_pairs, spacing=defspacing)
        j = j.ravel()
        j = j.tolist()
        for i in xrange(len(j)):
            matrix_to_list.append(j[i])
            originalJac.append(j[i])
    else:
        data = data.ravel()
        data = data.tolist()
        data = npr.choice(data, len(data))
        data.shape = shape
        print data.shape
        j, je, sf, sfe = calculate_average_Jacobian(model,fitopts, FRET_pairs=def_FRET_pairs, spacing=defspacing)
        j = j.ravel()
        j = j.tolist()
        for i in xrange(len(j)):
            matrix_to_list.append(j[i])
#    print len(c)
print 'Analysis Complete'
print 'Reverting list of data into a matrix'
list_to_matrix = np.asarray(matrix_to_list)
list_to_matrix.shape = (brs+1,-1)
print list_to_matrix.shape

print 'Creating directory for  data in output_files dir'
os.chdir('/dascratch/mjh9/output_files')
print os.getcwd()
from time import strftime
date = strftime('%Y_%b_%d')
print date
os.mkdir('bootJac_'+args.tempfile+date)
os.chdir('bootJac_'+args.tempfile+date)

print 'Calculating and appending the standard deviations of the data'
#print c.shape
stddevs = []
for i in xrange(brs+1):
    contact_i = [list_to_matrix[i,:]]
    stddevs.append(std(contact_i))
stdfile = open('std_T'+args.tempfile+str(brs)+'_resamples.txt', 'w+')
# list of standard deviations for each conctact, should be (num_bins*977) long
for i in xrange(len(stddevs)):
    stdfile.write(str(stddevs[i]))
    stdfile.write(' ')

print 'Calculating and appending the average values of the Jacobians'
avgJac = []
new_matrix = list_to_matrix.ravel()
new_matrix.shape = (brs+1, num_bins, 977)
print new_matrix.shape
for i in xrange(977):
    for j in xrange(num_bins):
        avg = mean(new_matrix[:,j,i])
        avgJac.append(avg)
avgfile = open('avgJac_T'+args.tempfile+str(brs)+'_resamples.txt', 'w+')
# this file will need to be converted into an array of dimensions (num_bins,977)
for i in xrange(len(avgJac)):
    avgfile.write(str(avgJac[i]))
    avgfile.write(' ')

print 'Arranging the Jacobians into a list'
final_list = new_matrix.tolist()
Jacfile = open('Jac_T'+args.tempfile+str(brs)+'_resamples.txt', 'w+')
# convert data into an array of dimensinos (brs+1, num_bins, 977)
for i in xrange(len(final_list)):
    Jacfile.write(str(final_list[i]))
    Jacfile.write(' ')

print 'Savings the original Jacobian'
orfile = open('originalJac_T'+args.tempfile+str(brs)+'_resamples.txt', 'w+')
for i in xrange(len(originalJac)):
    orfile.write(str(originalJac[i]))
    orfile.write(' ')

print 'COMPELETE'
exit()
