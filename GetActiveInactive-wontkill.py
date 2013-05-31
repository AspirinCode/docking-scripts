from fitting import *
import random
import optparse
import pylab
import numpy
import os
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from numpy import linalg

#bi path score samples
#[10,   1,     2,   4,   8,   3,                5]
#[2.0, 3.0,   5.0, 6.0, 7.0, 8.0,              12.0]

#car path score samples
#[39, 15,  14,     19,   5,   20,   3,   4]
#[2.0, 3.0,4.0,     6.0, 7.0, 8.0,   9.0, 11.0]

#apo path score samples
#[23,  108, 150, 71,  21,  7,        21,   20]
#[1.0, 2.0, 3.0, 4.0, 5.0, 6.0,      8.0,  9.0]

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)
def bootstrap(  b, n ):
    #Randomly divide n data points into b blocks.
    s = [ random.randint( 0, n-1 ) for t in xrange(0, b) ]
    return numpy.array(s)

def bar_progress(data, aucname):
    systems=['bi', 'car', 'apo']
    colors=['red', 'blue', 'grey' ]
    names=['inactive', 'inter1', 'inter2', 'active']
    index=0
    fig=pylab.figure(figsize=(8,6))
    #gs=gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1=fig.add_subplot(311, aspect=5)
    ax2=fig.add_subplot(312, aspect=5)
    ax3=fig.add_subplot(313, aspect=5)
    ax=[ax1, ax2, ax3]
    width=0.8
    for (n, sys) in enumerate(systems):
        avgs=[]
        errs=[]
        for key in names:
            avgs.append(numpy.mean(data[sys][key])) 
            errs.append(numpy.std(data[sys][key]))
        print avgs, errs
        ax[n].bar(left=range(0, len(names)), height=avgs,  alpha=0.5, width=width, yerr=errs,  ecolor='k')
    if aucname=='types':
        limits=(0.4,1.0)
    else:
        limits=(0.6,1.0)
    ax1.set_ylim(limits[0], limits[1])
    ax2.set_ylim(limits[0], limits[1])
    ax3.set_ylim(limits[0], limits[1])
    pylab.setp( ax1.get_xticklabels(), visible=False)
    pylab.setp( ax2.get_xticklabels(), visible=False)
    #ax1.set_ylabel('AUCs')
    ax1.yaxis.set_ticks(numpy.arange(0.6,1.05,0.1))
    ax2.yaxis.set_ticks(numpy.arange(0.6,1.05,0.1))
    ax3.yaxis.set_ticks(numpy.arange(0.6,1.05,0.1))
    ax3.xaxis.set_ticks(numpy.arange(0,4))
    ax3.xaxis.set_ticklabels(['Inactive', ' ', ' ', 'Active'])
    pylab.show()

def modify_scores(allscores):
    combine_inds=dict()
    combine_inds['bi']=[[2], [3, 5,6],[7,8], [12,]] 
    combine_inds['car']=[[2,], [3,4,6],[7,8], [9,11]] 
    combine_inds['apo']=[[1,2], [3,4,5,6],[8,], [9]] 
    newscores=dict()
    for sys in ['bi', 'car', 'apo']:
        newscores[sys]=dict()
        names=['inactive', 'inter1', 'inter2', 'active']
        print sys
        for (n,int) in enumerate(combine_inds[sys]):
            print names[n]
            if names[n] not in newscores[sys].keys():
                newscores[sys][names[n]]=[]
            for i in int:
                for val in allscores[sys][i]:
                    newscores[sys][names[n]].append(val)
    return newscores

def wrapper(func, args):
    c, p=func(*args)
    return c,p

def plot_fits(scores, values, R, pval, format, aucname):
    if pval< 0.0001:
        pval=0.0001
    pylab.figure()
    pylab.plot(scores, values, format)
    (ar,br)=polyfit(scores, values, 1)
    xr=polyval([ar,br], scores)
    pylab.plot(scores,xr,'%s-' % format[0], label='R=%s, pval=%s' %
            (round(R,2), round(pval,4)))
    #a, b, sa, sb, rchi2, dof=linear_fit(numpy.array(sorted(histo.keys())), avgs, stds)
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
    pylab.ylabel('%s Aucs' % (auclabels[aucname]))
    pylab.savefig('%s_%saucs.png' % (sys, aucname), dpi=300)

def get_fried(data):
    chis=[]
    pvals=[]
    cutoff=20
    for x in range(0,10):
        subsamples=dict()
        for key in data.keys():
            indices=bootstrap( cutoff, len(data[key]))
            subsamples[key]=data[key][indices]
        chi2, pval=wrapper(stats.friedmanchisquare, [subsamples[key] for key in subsamples.keys()])
        chis.append(chi2)
        pvals.append(pval)
    print "average chi2 %s+/%s:" % (numpy.mean(chis), numpy.std(chis))
    print "average pval %s+/-%s:" % (numpy.mean(pvals), numpy.std(pvals))
    return numpy.mean(chis), numpy.std(chis), numpy.mean(pvals), numpy.std(pvals)


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

def get_scores():
    max_score=0
    op_scores=dict()
    op_scores['h36']=dict()
    op_scores['conn']=dict()
    op_scores['npxxy']=dict()
    op_scores['bulge']=dict()
    vals=numpy.arange(8, 13,1)
    scores=range(0,len(vals))
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['h36'][v]=s
    vals=numpy.arange(0.5, 3.0, 0.5)
    scores=range(0,len(vals))[::-1]
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['conn'][v]=s
    vals=numpy.arange(0.5, 2.0, 0.5)
    scores=range(0,len(vals))
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['npxxy'][v]=s
    vals=numpy.arange(0.5, 2.5, 0.5)
    scores=range(0,len(vals))[::-1]
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['bulge'][v]=s
    return max_score, op_scores

def main(aucname):
    systems=['bi', 'car', 'apo']
    labels=['Agonist-bound', 'Inv. Agonist-bound', 'Apo']
    auclabels=dict()
    auclabels['types']='Discrimination'
    auclabels['antagonist']='Antagonist'
    auclabels['agonist']='Agonist'

    formats=['ro', 'bo', 'ko']
    allvalues=dict()
    allscores=dict()
    for (sys, label, format) in zip(systems, labels, formats):
        allscores[sys]=dict()
        dir='./%s/structural-data/' % sys
        pops=numpy.loadtxt('%s/Populations.dat' % dir)
        types=['path', 'new-path']
        type='path'
        states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(0,), dtype=int)
        pops=pops[states]
        if sys=='bi':
            cutoff=min(pops)
        else:
            frames=numpy.where(pops>= cutoff)[0]
            pops=pops[frames]
            states=states[frames]
        aucs=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(1,))
        paths=open('%s/%s_paths.txt' % (sys, sys))
        map=numpy.loadtxt('%s/Mapping.dat' % dir)
        pop=numpy.loadtxt('%s/Populations.dat' % dir)
        ops=dict()
        all_aucs=dict()
        max_score, op_scores=get_scores()
        print "max score is ", max_score
        all_active=[]
        all_inactive=[]
        for (i, path) in enumerate(paths.readlines()):
            path=numpy.array(path.split())
            score=numpy.zeros(len(path))
            h36=numpy.loadtxt('%s/%s%s.h36.dat' % (dir, type, i))
            bulge=numpy.loadtxt('%s/active-h5buldge-%s%s.dat' % (dir, type, i))
            in_npxxy=numpy.loadtxt('%s/inactive-npxxy-%s%s.dat' % (dir, type, i))
            ac_npxxy=numpy.loadtxt('%s/active-npxxy-%s%s.dat' % (dir, type, i))
            in_conn=numpy.loadtxt('%s/inactive-conn-%s%s.dat' % (dir, type, i))
            ac_conn=numpy.loadtxt('%s/active-conn-%s%s.dat' % (dir, type, i))
            score=op_path(h36, 'h36', op_scores, score)
            score=op_path(ac_conn, 'conn', op_scores, score)
            score=op_path(in_npxxy, 'npxxy', op_scores, score)
            score=op_path(bulge, 'bulge', op_scores, score)
            ops[i]=[]
            all_aucs[i]=[]
            for (m, state) in enumerate(path):
                location=numpy.where(states==int(state))[0]
                if location.size:
                    all_aucs[i].append(aucs[location][0])
                    ops[i].append(score[m])
                    if score[m] not in allscores[sys].keys():
                        allscores[sys][score[m]]=[]
                        allscores[sys][score[m]].append(aucs[location][0])
                    else:
                        allscores[sys][score[m]].append(aucs[location][0])
        print "building score hist"
        histo=dict()
        for path in sorted(ops.keys()):
            for (n, score) in enumerate(ops[path]):
                if score not in histo.keys():
                    histo[score]=[]
                    histo[score].append(all_aucs[path][n])
                else:
                    histo[score].append(all_aucs[path][n])

        values=[]
        scores=[]
        for s in sorted(histo.keys()):
            for i in histo[s]:
                values.append(i)
                scores.append(s)
        values=numpy.array(values)
        scores=numpy.array(scores)
        slope, intercept, R, pval, std_err = stats.linregress(scores,  values)
        print sys, R, pval
        #plot_fits(scores, values, R, pval, format, aucname)
        allvalues[sys]=values
    chi2, chi2_s, pval, pval_s=get_fried(allvalues)
    print "on pathway score and AUCS"
    #for sys in systems:
    #    for s in range(1,13):
    #      not in allscores[sys].keys():
    #            allscores[sys][s]=0
    newscores=modify_scores(allscores)
    bar_progress(newscores, aucname)
    #chi2, chi2_s, pval, pval_s=get_fried(newscores)


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-n', '--aucname', dest='aucname',
                      help='auc type name')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(aucname=options.aucname)


