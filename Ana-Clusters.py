from DockingStats import *
import glob
import random
import optparse
import pylab
import numpy
import os
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from numpy import linalg


def main( dir, keynumber): 
    keys=['h36', 'conn', 'npxxy', 'bulge']
    keys=keys[:int(keynumber)]
    print "Using OP's :", keys
    files=glob.glob('%s/rocs-bw_*-*combo*g*dat' % dir)
    #files=numpy.loadtxt('cluster-results/agonist/top0.5_populated_clusters.txt',
    #        dtype=str)
    names=[]
    counts=[]
    ligands=['bi',  'car', 'apo']
    topdata=dict()
    #path scores for reference
    max_score, op_scores=get_ref_scores(keys)
    for file in files:
        name=file.split('bw_')[1].split('_combo')[0]
        sys=file.split('bw_')[1].split('-')[0]
        state=file.split('bw_')[1].split('-')[1]
        if sys not in topdata.keys():
            topdata[sys]=[]
        data=numpy.loadtxt(file, usecols=(0,), dtype=str)
        if data.size > 23:
            topdata[sys].append(state)
            names.append(name)
            counts.append(data.size)
            subdata=dict()
            for i in range(0, data.size):
                if data.size==1:
                    x=str(data)
                else:
                    x=str(data[i])
                sys=x.split('-')[0]
                state=int(x.split('-')[1])
                if sys not in subdata.keys():
                    subdata[sys]=[]
                subdata[sys].append(state)
                pop=numpy.loadtxt('%s/cent-structural-data/Populations.dat' % sys)[state]
            #for sys in subdata.keys():
            #    SD=SystemDock(sys, subdata[sys])
    names=numpy.array(names)
    counts=numpy.array(counts)
    frames=numpy.argsort(counts)
    for file in names[frames][::-1]: # sorted from most to least members
        location=numpy.where(names==file)[0]
        sys=file.split('-')[0]
        state=int(file.split('-')[1])
        SD=SystemDock(sys, numpy.array([state,]))
        #SD=SystemDock(sys, topdata[sys])
        pathfile=open('%s/%s_paths.txt' % (sys, sys))
        SD.get_path_scores('%s/cent-structural-data/' % sys, pathfile, keys, op_scores)
        scores=SD.pathscores
        print '%s-%s' % (sys, state), "progress scores: ", scores, "%s members" % counts[location]


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
            help='directory for gen dat files')
    parser.add_option('-k', '--keynumber', dest='keynumber',
            help='# OPs (1- bulge, 2 - conn, 3 - npxxy, 4 - h36')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, keynumber=options.keynumber)
