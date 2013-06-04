from fitting import *
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
def bootstrap(  b, n ):
    #Randomly divide n data points into b blocks.
    s = [ random.randint( 0, n-1 ) for t in xrange(0, b) ]
    return numpy.array(s)

def get_aggregate(scoredata, aggdata, keys, limit=False):
    if limit==True:
        print "only getting aggregate ligand data"
        systems=['bi', 'car']
    else:
        systems=['bi', 'car', 'apo']
    for s in systems:
        for key in keys:
            try: 
                test=int(key)
                if test>=9:
                    key=9
            except ValueError:
                pass
            if key in scoredata[s].keys():
                for val in scoredata[s][key]:
                    if key not in aggdata.keys():
                        aggdata[key]=[]
                    aggdata[key].append(val)
            else:
                pass
    return aggdata, aggdata.keys()



def bar_progress(data):
    auclabels=dict()
    auclabels['types']='      Discrimination'
    auclabels['antagonist']='Antagonist'
    auclabels['agonist']='Agonist'

    systems=['bi', 'car', 'apo']
    labels=['Agonist-bound', 'Inv. Agnonist bound', 'Apo']
    colors=['red', 'blue', 'grey' ]
    names=['inactive', 'inter', 'active']
    index=0
    fig=pylab.figure(figsize=(8,6))
    #gs=gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1=fig.add_subplot(321, aspect=7)
    ax2=fig.add_subplot(322, aspect=7)
    ax3=fig.add_subplot(323, aspect=7)
    ax4=fig.add_subplot(324, aspect=7)
    ax5=fig.add_subplot(325, aspect=7)
    ax6=fig.add_subplot(326, aspect=7)
    ax=[ax1, ax3, ax5, ax2, ax4, ax6]
    width=0.8
    count=0
    for aucname in ['agonist', 'antagonist']:
        limits=(0.6,1.0)
        for (n, sys) in enumerate(systems):
            avgs=[]
            errs=[]
            for key in names:
                avgs.append(numpy.mean(data[aucname][sys][key])) 
                errs.append(numpy.std(data[aucname][sys][key]))
            if aucname=='agonist':
                ax[count].bar(left=range(-3,-2), height=[0,], color=colors[n],
                        label=labels[n])
                lg=ax[count].legend(loc=9)
                lg.draw_frame(False)

            ax[count].bar(left=range(0, len(names)), height=avgs,  
                    width=width, yerr=errs,  color=colors[n], ecolor='k')
            ax[count].set_ylim(limits[0], limits[1])
            ax[count].yaxis.set_ticks(numpy.arange(limits[0],limits[1]+0.05,0.1))
            ax[count].xaxis.set_ticks(numpy.arange(0,4))
            ax[count].xaxis.set_ticklabels([' ', 'Inactive',  ' ', 'Active'])
            count+=1
    pylab.setp( ax1.get_xticklabels(), visible=False)
    pylab.setp( ax2.get_xticklabels(), visible=False)
    pylab.setp( ax3.get_xticklabels(), visible=False)
    pylab.setp( ax4.get_xticklabels(), visible=False)
    pylab.text(0.5, 1.15, 'Agonist vs. Decoys',
                     horizontalalignment='center',
                              fontsize=16, transform = ax1.transAxes)
    pylab.text(0.5, 1.15, 'Antagonist vs. Decoys',
                     horizontalalignment='center',
                              fontsize=16, transform = ax2.transAxes)
    ax3.set_ylabel('AUCs')
    ax5.set_xlabel('Path Progress')
    ax6.set_xlabel('Path Progress')
    fig.subplots_adjust(wspace=0.001)
    pylab.savefig('ligand_path_aucs.png', dpi=300)
    ###########
    aucname='types'
    limits=(0.4,1.0)
    fig=pylab.figure()
    ax1=fig.add_subplot(311, aspect=5)
    ax2=fig.add_subplot(312, aspect=5)
    ax3=fig.add_subplot(313, aspect=5)
    count=0
    ax=[ax1, ax2, ax3]
    for (n, sys) in enumerate(systems):
        avgs=[]
        errs=[]
        for key in names:
            avgs.append(numpy.mean(data[aucname][sys][key])) 
            errs.append(numpy.std(data[aucname][sys][key]))
        if sys=='bi':
            ax[count].plot(range(0,5), [0.5]*5, 'k--', label='Random')
            lg=ax1.legend(loc=9)
            lg.draw_frame(False)
        else:
            ax[count].plot(range(0,5), [0.5]*5, 'k--', label='Random')
        ax[count].bar(left=range(0, len(names)), height=avgs,  
                width=width, yerr=errs, color=colors[n], ecolor='k')
        ax[count].set_ylim(limits[0], limits[1])
        ax[count].yaxis.set_ticks(numpy.arange(limits[0],limits[1]+0.05,0.1))
        ax[count].xaxis.set_ticks(numpy.arange(0,3))
        ax[count].xaxis.set_ticklabels([' ', 'Inactive',  ' ', 'Active'])
        count+=1
    pylab.text(0.5, 1.15, 'Agonist vs. Inv. Ag.',
                     horizontalalignment='center',
                              fontsize=16, transform = ax1.transAxes)
    pylab.setp( ax1.get_xticklabels(), visible=False)
    pylab.setp( ax2.get_xticklabels(), visible=False)
    ax3.set_xlabel('Path Progress')
    ax2.set_ylabel('AUCs')
    #ax1.set_ylabel('AUCs')
    #pylab.savefig('types_path_aucs.png', dpi=300)
    pylab.show()

def modify_scores(all, new):
    combine_inds=dict()
    combine_inds['bi']=[[2], [3, 5],[6,7], [8,12,]] 
    combine_inds['car']=[[2,], [3,4],[6,7], [8, 9,11]] 
    combine_inds['apo']=[[1,2], [3,4],[5,6], [8,9]] 
    #combine_inds['bi']=[[2,3], [5, 6,7], [8,12]] 
    #combine_inds['car']=[[2,3 ], [4, 6, 7], [8, 9,11]] 
    #combine_inds['apo']=[[1,2,3 ], [4, 5,6], [8,9]] 
    if len(combine_inds['bi'])==3:
        names=['inactive', 'inter', 'active']
    elif len(combine_inds['bi'])==4:
        names=['inactive', 'inter1', 'inter2',  'active']
    new=dict()
    for sys in ['bi', 'car', 'apo']:
        new[sys]=dict()
        for (n,int) in enumerate(combine_inds[sys]):
            if names[n] not in new[sys].keys():
                new[sys][names[n]]=[]
            for i in int:
                for val in all[sys][i]:
                    new[sys][names[n]].append(val)
        for key in names:
            new[sys][key]=numpy.array(new[sys][key])
    return new

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

def get_fried(data, keys):
    chis=[]
    pvals=[]
    cutoff=20
    #if len(keys) < 3:
    #    chi2, pval=wrapper(stats.chisquare, [data[key] for key in keys])
    k_chi2, k_pval=wrapper(stats.mstats.kruskalwallis, [data[key] for key in keys])
    print "kruskal chi2 %s p val %s :" % (k_chi2, k_pval)
    return numpy.mean(chis), numpy.std(chis), numpy.mean(pvals), numpy.std(pvals), k_chi2, k_pval


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

def get_scores(opnum):
    max_score=0
    op_scores=dict()
    op_scores['h36']=
    op_scores['conn']=dict()
    op_scores['npxxy']=dict()
    op_scores['bulge']=dict()

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

def main():
    labels=['Agonist-bound', 'Inv. Agonist-bound', 'Apo']
    aucnames=['agonist', 'antagonist', 'types']
    auclabels=dict()
    auclabels['types']='      Discrimination'
    auclabels['antagonist']='Antagonist'
    auclabels['agonist']='Agonist'
    newscores=dict()
    allscores=dict()
    allvalues=dict()
    for aucname in aucnames:
        allscores[aucname]=dict()
        allvalues[aucname]=dict()
        newscores[aucname]=dict()
        systems=['bi', 'car', 'apo']
        labels=['Agonist-bound', 'Inv. Agonist-bound', 'Apo']
        auclabels=dict()
        auclabels['types']='Discrimination'
        auclabels['antagonist']='Antagonist'
        auclabels['agonist']='Agonist'
        formats=['ro', 'bo', 'ko']
        for (sys, label, format) in zip(systems, labels, formats):
            allscores[aucname][sys]=dict()
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
            statetracker=[]
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
                    if state not in statetracker:
                        statetracker.append(state)
                        location=numpy.where(states==int(state))[0]
                        if location.size:
                            all_aucs[i].append(aucs[location][0])
                            ops[i].append(score[m])
                            if score[m] not in allscores[aucname][sys].keys():
                                allscores[aucname][sys][score[m]]=[]
                                allscores[aucname][sys][score[m]].append(aucs[location][0])
                            else:
                                allscores[aucname][sys][score[m]].append(aucs[location][0])
                    else:
                        pass
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
            print sys, "correlation path scores, auc values: ",  R, pval
            #plot_fits(scores, values, R, pval, format, aucname)
            allvalues[aucname][sys]=values
        print "Chi2 Test for All Values %s" % aucname
        if aucname=='types':
            print "using only ligand data"
            systems=['bi', 'car']
        else:
            systems=['bi', 'car', 'apo']
        chi2, chi2_s, pval, pval_s, k_chi2, k_pval=get_fried(allvalues[aucname], systems)
        print "on pathway score and AUCS"
        #for sys in systems:
        #    for s in range(1,13):
        #      not in allscores[aucname][sys].keys():
        #            allscores[aucname][sys][s]=0
        newscores[aucname]=modify_scores(allscores[aucname], newscores[aucname])
    import pdb
    pdb.set_trace()
    names=['inactive', 'inter', 'active']
    #bar_progress(newscores, auclabels)
    allprogress=dict()
    limit=False
    for aucname in aucnames:
        for s in systems:
            print "Chi2 Test for %s AUCs and %s Path States" % (aucname, s)
            chi2, chi2_s, pval, pval_s, k_chi2, k_pval=get_fried(newscores[aucname][s], names)
        ohandle=open('%s_reduced_allvalues.dat' % aucname, 'w')
        allprogress=dict()
        if aucname=='types':
            limit=True
        else:
            limit=False
        allprogress, modscores=get_aggregate(newscores[aucname], allprogress, names, limit)
        print "Chi2 Test for %s AUCs and All Path Reduced States" % (aucname)
        for key in sorted(allprogress.keys()):
            for value in allprogress[key]:
                ohandle.write('%s\t%s\n' % (key, value)) 
        chi2, chi2_s, pval, pval_s, k_chi, k_pval=get_fried(allprogress, names)

        ohandle=open('%s_allvalues.dat' % aucname, 'w')
        allprogress=dict()
        allprogress, modscores=get_aggregate(allscores[aucname], allprogress, range(1,13), limit)
        print "Chi2 Test for %s AUCs and All Path States" % (aucname)
        for key in sorted(allprogress.keys()):
            for value in allprogress[key]:
                ohandle.write('%s\t%s\n' % (key, value)) 
        chi2, chi2_s, pval, pval_s, k_chi, k_pval=get_fried(allprogress, modscores)

if __name__ == "__main__":
    main()


