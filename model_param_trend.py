import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import argparse
import os

parser = argparse.ArgumentParser(description = "enter the path, file name; will then create graph of epsilon's progression throughout each iteration")
parser.add_argument('path', help = 'data for each temperature can be found through this path; can also simply lead straight to the desired file')
parser.add_argument('file', help = 'name of the final file required from the data in each temp dir, include / first')
args = parser.parse_args()

fig = plt.figure()
ax = fig.add_subplot(111)

###################################PRESET VARIABLES#######################################
action = 'trend'
#for trend
plot = 'all'
filter_switch = 'y'
min_change = 0.5
#for limit
all_or_nothing = ''
#for range
a = ''
b = ''
#for one
n = ''
##########################################################################################

 # action = raw_input('Which code would you like to run? trend/limit ')
action = 'trend'
if action == 'trend':
#TREND
 #    plot = raw_input('Would you like to plot all contacts, just one, or a range? answer "all", "range", or "one" ')
    if plot == 'one':
 #        n = input('Which contact would you like to see? ')
        for x in xrange(n,n+1):
            #this is the list to be filled with the ten specific values and then later graphed
            variable_x=[1]
            for i in xrange(args.tmin, args.tmax, args.tstep):
                #this will keep the code from adding more than 10 values to the list and from trying to open files that don't exist
                #opening the files
                model_params_i = np.loadtxt(args.path+str(i)+args.file)
                #turning the data set into an array
                np.asarray(model_params_i)
                #appends one value from each file into the previously created list
                variable_x.append(model_params_i[x])
            #labels for the graph
#            plt.xlabel('Iterations')
#            plt.ylabel('Model Parameters')
            #this next part controls whether or not a filter will be applied to the results
 #            filter_switch=raw_input('Would you like to filter out changes in parameter if it is below a certain number? y/n ')
            if filter_switch == 'y':
 #                min_change = float(input('Enter the minimum change: '))
                delta = []
                try:
                    for z in xrange(10):
                        delta.append(abs(1-variable_x[z]))
                        chg = max(delta)
                        #this part is the actual filter
                    if chg < min_change:
                        variable_x = []
                except IndexError:
                    pass
                #time to plot!
                ax.plot(variable_x, 'b')
#                plt.title('Filtered '+str(min_change)+', Contact '+str(x))
            elif filter_switch == 'n':
                ax.plot(variable_x, 'g')
#                plt.title('Unfiltered, Contact '+str(x))
            else:
                print 'Please enter y or n'
                exit()
            fig.savesig('contact_'+str(n)+'.png')
    elif plot == 'all':
 #        filter_switch=raw_input('Would you like to filter out changes in parameter if it is below a certain number? y/n ')
 #        if filter_switch == 'y':
 #            min_change = float(input('Enter the minimum change: '))
 #        elif filter_switch == 'n':
 #            pass
 #        else:
 #            print 'Please enter y or n'
        for x in xrange(977):
            variable_x=[1]
            for i in xrange(10):
                model_params_i = np.loadtxt(args.path+str(i)+args.file)
                np.asarray(model_params_i)
                variable_x.append(model_params_i[x])
                #print variable_x
#            plt.xlabel('Iterations')
#            plt.ylabel('Model Parameters')
            if filter_switch=='y':
                delta = []
                try:
                    for z in xrange(10):
                        delta.append(abs(1-variable_x[z]))
                        chg = max(delta)
                    if chg < min_change:
                        variable_x = []
                except IndexError:
                    pass
                ax.plot(variable_x, 'b')
#                plt.title('Filtered '+str(min_change)+', All Contacts')
            else:
                ax.plot(variable_x, 'g')
#                plt.title('Unfiltered, All Contacts')
        fig.savefig('all_contacts.png')
    elif plot == 'range':
 #        a = input('At which contact number would you like to begin the analysis? ')
 #        b = input('At which contact number would you like to end the analysis? (exclusive) ')
 #        filter_switch = raw_input('Filter? y/n ')
 #        if filter_switch == 'y':
 #            min_change=float(input('Enter the minimum change: '))
 #        elif filter_switch == 'n':
 #            pass
 #        else:
 #            print 'Please enter y or n'
        for x in xrange(a,b):
            variable_x=[1]
            for i in xrange(10):
                model_params_i = np.loadtxt(args.path+str(i)+args.file)
                variable_x.append(model_params_i[x])
#            plt.xlabel('Iterations')
#            plt.ylabel('Model Parameters')
            if filter_switch=='y':
                delta = []
                try:
                    for z in xrange(10):
                        delta.append(abs(1-variable_x[z]))
                        chg = max(delta)
                    if chg < min_change:
                        variable_x = []
                except IndexError:
                    pass
                ax.plot(variable_x, 'b')
#                plt.title('Filtered '+str(min_change)+', Contacts '+str(a)+' to '+str(b))
            else:
                ax.plot(variable_x, 'g')
#                plt.title('Unfiltered, Contacts '+str(a)+' to '+str(b))
        fig.savefig('contacts_'+str(a)+'-'+str(b)+'.png')
    else:
        print 'Please enter a valid answer'
elif action == 'limit':
#LIMIT
 #    all_or_nothing = raw_input('Which contacts would you like to analyze? answer "one", "all", "none", or "range" ')
#    plt.xlabel('Iteration')
#    plt.ylabel('Model Parameter')
    if all_or_nothing == 'all':
        contact_num = []
        for t in xrange(977):
            trace_t = [1]
            for i in xrange(10):
                data_i = np.loadtxt(args.path+str(i)+args.file)
                #print data_i[t]
                trace_t.append(data_i[t])
                if len(trace_t) == 10:
                    x = min(trace_t)
                    if x > 0.01:
                        trace_t = []
                    else:
                        ax.plot(trace_t)
                        contact_num.append(str(t))
                        #print t
#            plt.title('Contacts Reaching 0.01')
        fig.savefig('all_contacts_reaching_0.01.png')
        #print contact_num
        print len(contact_num)
    elif all_or_nothing == 'one':
 #        n = input('Which contact would you like to analyze? ')
        for t in xrange(n,n+1):
            trace_t = [1]
            for i in xrange(10):
                data_i = np.loadtxt(args.path+str(i)+args.file)
                #print data_i[t]
                trace_t.append(data_i[t])
                if len(trace_t) == 10:
                    x = min(trace_t)
                    if x > 0.01:
                        trace_t = []
                    else:
                        ax.plot(trace_t)
                        #print t
#            plt.title('Contacts Reaching 0.01')
        fig.savefig('contact_'+str(n)+'_reaching_0.01.png')
    elif all_or_nothing == 'range':
        contact_num = []
 #        a = input('At which contact would you like to begin analysis? ')
 #        b = input('At which contact would you like to end analysis? ')
        for t in xrange(a,b):
            trace_t = [1]
            for i in xrange(10):
                it_num = i
                data_i = np.loadtxt(args.file+str(it_num))
                #print data_i[t]
                trace_t.append(data_i[t])
                if len(trace_t) == 10:
                    x = min(trace_t)
                    if x > 0.01:
                        trace_t = []
                    else:
                        ax.plot(trace_t)
                        contact_num.append(str(t))
                        #print t
#            plt.title('Contacts Reaching 0.01')
        fig.savefig('contacts_'+str(a)+'-'+str(b)+'_reaching_0.01.png')
        #print contact_num
        print len(contact_num)
    elif all_or_nothing == 'none':
        pass
    #Part 2
    pairwise = np.loadtxt(args.file)
    param = np.asarray(contact_num)
    param = param.astype(float)
    contacts = zeros((51,2),dtype=float)
    for i in xrange(51):
        num = param[i]
        #print num
        for t in xrange(977):
            residues = pairwise[t,2]
            #print residues
            if residues == num:
                print pairwise[t,:2]
                contacts[i,:2] = pairwise[t,:2]

exit()
