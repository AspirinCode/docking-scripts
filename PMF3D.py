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

def wrapper(func, args):
    c, p=func(*args)
    return c,p

def get_chi2(data):
    if len(data) < 3:
        print "chi2 test for two samples"
        chi2, pval=wrapper(stats.chisquare, data)
    else:
        chi2, pval=wrapper(stats.mstats.kruskalwallis, data)
    print "kruskal chi2 %s p val %s :" % (chi2, pval)
    return chi2, pval

def get_ref_scores(keys):
    max_score=0
    ref=dict()
    vals=numpy.arange(8, 13,1)
    ref['h36']=vals
    vals=numpy.arange(0.5, 3.0, 0.5)
    ref['conn']=vals
    vals=numpy.arange(8, 13,1)
    ref['npxxy']=vals
    vals=numpy.arange(0.5, 2.5, 0.5)
    ref['bulge']=vals
    op_scores=dict()
    for key in keys:
        op_scores[key]=dict()
        if key=='bulge' or key=='conn':
            scores=range(0,len(ref[key]))[::-1]
        else:
            scores=range(0,len(ref[key]))
        max_score+=max(scores)
        for (v,s) in zip(ref[key], scores):
            op_scores[key][v]=s
    return max_score, op_scores

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

def get_label(name):
    if name=='types':
        label='Agonist vs. Antag.'
    elif name=='antagonist':
        label='Antagonist vs. Decoy'
    elif name=='agonist':
        label='Agonist vs. Decoy'
    elif name=='bi':
        label='Agonist-bound'
    elif name=='car':
        label='Inv. Agonist bound'
    elif name=='apo':
        label='Apo'
    return label


def scatter_fits(sys, scores, values, R, pval, aucname, show=False):
    auclabel=get_label(aucname)
    if sys=='apo':
        format='ko'
        label=
    if sys=='bi':
        format='ro'
    if sys=='car':
        format='bo'
    if pval< 0.0001:
        pval=0.0001
    pylab.figure()
    pylab.plot(scores, values, format)
    (ar,br)=polyfit(scores, values, 1)
    xr=polyval([ar,br], scores)
    pylab.plot(scores,xr,'%s-' % format[0], label='R=%s, pval=%s' %
            (round(R,2), round(pval,4)))
    if aucname=='types':
        pylab.plot(range(0, 15), [0.5]*len(range(0,15)), 'k--', label='Random Disc.')
        pylab.ylim(0.3, 0.9)
    else:
        pylab.ylim(0.5, 1.0)
    pylab.xlim(0, 15)
    lg=pylab.legend()
    lg.draw_frame(False)
    pylab.title('%s States' % label)
    pylab.xticks(range(0, 15), [' ']*2+ ['inactive']+[' ']*(len(range(0,15))-6)+['active']+[' ']*2)
    pylab.xlabel('Pathway Progress')
    pylab.ylabel('%s Aucs' % (auclabel))
    pylab.savefig('%s_%saucs.png' % (sys, aucname), dpi=300)
    if show==True:
        pylab.show()

# Class for AUC Data
class SystemDock:
    def __init__(self, ligand, states, aucs, pops):
        self.ligand=ligand
        self.aucs = aucs
        self.states = states
        self.pops = pops
        self.pathscores=[]
        self.reducescores=[]
        self.mapscores=[]

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

    def get_path_scores(dir, pathfile, keys):
        statescores=numpy.zeros(len(SystemDock.states))
        for (i, path) in enumerate(pathfile.readlines()):
            path=numpy.array(path.split())
            frames=[]
            for state in path:
                if state in SystemDock.states:
                    frames.append(state)
            frames=numpy.array(frames)
            path_ops=dict()
            op_scores=dict()
            for key in keys:
                if key=='h36':
                    path_ops[key]=numpy.loadtxt('%s/path%s.%s.dat' % (dir, key, i))
                    path_ops[key]=path_ops[key][frames]
                elif key=='npxxy':
                    path_ops[key]=numpy.loadtxt('%s/inactive-%s-path%s.dat' % (dir, key, i))
                    path_ops[key]=path_ops[key][frames]
                elif key=='conn':
                    path_ops[key]=numpy.loadtxt('%s/active-%s-path%s.dat' % (dir, key, i))
                    path_ops[key]=path_ops[key][frames]
                elif key=='bulge':
                    path_ops[key]=numpy.loadtxt('%s/active-h5%s-path%s.dat' % (dir, key, i))
                    path_ops[key]=path_ops[key][frames]
            score=numpy.zeros(len(path[frames]))
            for key in keys:
               score=op_path(path_ops[key], key, op_scores, score)
            for (i, state) in enumerate(path):
                location=numpy.where(SystemDock.states==state)[0]
                statescores[location]=score[i]
       self.pathscores=statescores 

    def modify_scores(length):
        combine_inds=dict()
        if length==4:
            combine_inds['bi']=[[2], [3, 5],[6,7], [8,12,]]
            combine_inds['car']=[[2,], [3,4],[6,7], [8, 9,11]]
            combine_inds['apo']=[[1,2], [3,4],[5,6], [8,9]]
            names=['inactive', 'inter1', 'inter2',  'active']
        if length==3:
            combine_inds['bi']=[[2,3], [5, 6,7], [8,12]]
            combine_inds['car']=[[2,3 ], [4, 6, 7], [8, 9,11]]
            combine_inds['apo']=[[1,2,3 ], [4, 5,6], [8,9]]
            names=['inactive', 'inter', 'active']
        new=numpy.zeros(len(SystemDock.pathscores))
        map_score=dict()
        key=SystemDock.ligand
        for score in SystemDock.pathscores:
            index=numpy.where(SystemDock.pathscores==score)[0]
            for (n,indices) in enumerate(combine_inds[key]):
                if score in indices:
                    new[index]=n
                    map[n]=names[n]
        self.reducescores=new
        self.mapscores=map

