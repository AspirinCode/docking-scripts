from math import  *
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
def bootstrap( b, n ):
    """
    Randomly divide n data points into b blocks.
    """
    s = [ random.randint( 0, n-1 ) for t in xrange(0, b) ]
    return s


def ttest_finish(df,t):
    """Common code between all 3 t-test functions."""
    prob = stats.distributions.t.sf(numpy.abs(t), df) * 2  # use numpy.abs to get upper tail
    if t.ndim == 0:
        t = t[()]

    return t, prob

def ztest(teststat, nullvalue,  se, alpha =0.05, side=0):
    """
    Normal test of a sample statistic.
    Return:
       pvalue
       Z test statistic
       critical value(s)

    Arguments:
      teststat- test statistic
      se -    standard error of sample statistic
      alpha - significance level
      side  - -1 left-sided test
             0  double sided test
             +1 right-sided test
    """
    pnorm = stats.norm.cdf
    qnorm = stats.norm.ppf
    Ztest = (teststat-nullvalue)/se
    #print "Ztest=",  Ztest
    if side ==0:
       pvalue = pnorm(Ztest)
       if Ztest > 0.0:
           pvalue = 1.0 -pvalue
       pvalue *= 2.0
       zcrit1 = qnorm(alpha/2)
       zcrit2 = qnorm(1-alpha/2.0)
       return pvalue, Ztest, (zcrit1,zcrit2)
    elif side == -1:
       pvalue = pnorm(Ztest)
       zcrit = qnorm(alpha)
       return pvalue, Ztest, zcrit
    else:
       pvalue = 1- pnorm(Ztest)
       zcrit  = qnorm(1.0-alpha)
       return pvalue, Ztest, zcrit

def ttest_ind(a, b, axis=0, equal_var=True):
    a, b, axis = stats._support._chk2_asarray(a, b, axis)
    v1 = numpy.var(a, axis, ddof=1)
    v2 = numpy.var(b, axis, ddof=1)
    n1 = a.shape[axis]
    n2 = b.shape[axis]

    if (equal_var):
        df = n1 + n2 - 2
        svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / float(df)
        denom = numpy.sqrt(svar * (1.0 / n1 + 1.0 / n2))
    else:
        vn1 = v1 / n1
        vn2 = v2 / n2
        df = ((vn1 + vn2)**2) / ((vn1**2) / (n1 - 1) + (vn2**2) / (n2 - 1))

        # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
        # Hence it doesn't matter what df is as long as it's not NaN.
        df = numpy.where(numpy.isnan(df), 1, df)
        denom = numpy.sqrt(vn1 + vn2)

    d = numpy.mean(a, axis) - numpy.mean(b, axis)
    t = numpy.divide(d, denom)
    t, prob = ttest_finish(df, t)

    return t, prob

def twosampleproptest(x1,  n1,   x2,  n2,  alpha,  side):
    p1hat = float(x1)/n1
    p2hat = float(x2)/n2
    phat = float(x1 +x2)/(n1+n2)
    #print "p1hat, p2hat",  p1hat,  p2hat
    samplestat =  p1hat - p2hat
    se = sqrt(phat*(1.0-phat) * (1.0/n1 + 1.0/n2))
    print "se %s" % se
    print "%s-sided" % side
    return ztest(p2hat, p1hat,  se,  alpha,  side)

def check(data, xtal):
    count=0
    for x in data:
        if x > xtal:
            count+=1
    return count



def wrapper(func, args):
    c, p=func(*args)
    return c,p

def chi2_test(data):
    if len(data) < 3:
        print "chi2 test for two samples"
    chi2, pval=wrapper(stats.mstats.kruskalwallis, data)
    return chi2, pval

def get_ref_scores(keys):
    max_score=0
    ref=dict()
    vals=numpy.arange(8, 13,1)
    ref['h36']=vals
    vals=numpy.arange(0.5, 3.0, 0.5)
    ref['conn']=vals
    vals=numpy.arange(0.5, 3.0, 0.5)
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
        label=get_label(sys)
    if sys=='bi':
        format='ro'
        label=get_label(sys)
    if sys=='car':
        format='bo'
        label=get_label(sys)
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

    def op_path(self, op, name, op_scores, score):
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

    def get_path_scores(self, dir, pathfile, keys, op_scores):
        statescores=numpy.zeros(len(self.states))
        for (i, path) in enumerate(pathfile.readlines()):
            path=numpy.array(path.split())
            frames=[]
            for (n, state) in enumerate(path):
                if int(state) in self.states:
                    frames.append((n-1))
            frames=numpy.array(frames)
            path_ops=dict()
            for key in keys:
                if key=='h36':
                    path_ops[key]=numpy.loadtxt('%s/path%s.%s.dat' % (dir, i, key))
                    path_ops[key]=path_ops[key][frames]
                elif key=='npxxy':
                    path_ops[key]=numpy.loadtxt('%s/inactive-%s-path%s.dat' % (dir, key, i))
                    path_ops[key]=path_ops[key][frames]
                elif key=='conn':
                    path_ops[key]=numpy.loadtxt('%s/active-%s-path%s.dat' % (dir, key, i))
                    path_ops[key]=path_ops[key][frames]
                elif key=='bulge':
                    path_ops[key]=numpy.loadtxt('%s/active-h5buldge-path%s.dat' % (dir, i))
                    path_ops[key]=path_ops[key][frames]
            score=numpy.zeros(len(path[frames]))
            for key in keys:
               score=op_path(path_ops[key], key, op_scores, score)
            for (i, state) in enumerate(path[frames]):
                location=numpy.where(self.states==int(state))[0]
                if statescores[location]!=0:
                    if statescores[location]!=score[i]:
                        print "problem wiht scoring"
                        import pdb
                        pdb.set_trace()
                    else:
                        statescores[location]=score[i]
                else:
                    statescores[location]=score[i]
        self.pathscores=statescores 

    def modify_scores(self, length):
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
        new=numpy.zeros(len(self.pathscores))
        map_score=dict()
        key=self.ligand
        for score in self.pathscores:
            index=numpy.where(self.pathscores==score)[0]
            for (n,indices) in enumerate(combine_inds[key]):
                if score in indices:
                    new[index]=n
                    map[n]=names[n]
        self.reducescores=new
        self.mapscores=map


