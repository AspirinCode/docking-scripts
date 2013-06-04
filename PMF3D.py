import random
import optparse
import pylab
import numpy
import os
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from numpy import linalg


font = {'family' : 'sans-serif',
                'sans-serif':['Helvetica'],
                                        'size'   : 14}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 10,
                  'legend.linewidth': 2}
pylab.rcParams.update(params)


# Helper Functions
        GDfree[i,j,k]=max
    return GDfree

# Class for AUC Data
class SystemDock:
    def __init__(self, ligand, states, aucs, pops):
        self.ligand=ligand
        self.aucs = aucs
        self.states = states
        self.pops = pops

    def op_path(op, name, op_scores, score):
        test=sorted(op_scores[name].keys())
        for (n, x) in enumerate(op):
            for (m, i) in enumerate(test):
                if x >= test[m]:
                    if (m+1) >= len(test):
                        score[n]+=op_scores[name][i]
                    elif x < test[m+1]:
                        score[n]+=op_scores[name][i]
                        break
        return score

