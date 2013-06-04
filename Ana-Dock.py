from fitting import *
from DockingStats import *
import random
import optparse
import pylab
import numpy
import os
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from numpy import linalg


def main():
    labels=['Agonist-bound', 'Inv. Agonist-bound', 'Apo']
    aucnames=['agonist', 'antagonist', 'types']
    newscores=dict()
    allscores=dict()
    allvalues=dict()

    systems=dict()
    for aucname in aucnames:
        allscores[aucname]=dict()
        allvalues[aucname]=dict()
        newscores[aucname]=dict()
        ligands=['bi', 'car', 'apo']
        for sys in ligands:
            states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(0,), dtype=int)
            dir='./%s/structural-data/' % sys
            pops=numpy.loadtxt('%s/Populations.dat' % dir)
            pops=pops[states]
            if sys=='bi':
                cutoff=min(pops)
            else:
                frames=numpy.where(pops>= cutoff)[0]
                pops=pops[frames]
                states=states[frames]
            aucs=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(1,))
            system=SystemDock(sys, states, aucs, pops):
            systems.append(system)
            paths=open('%s/%s_paths.txt' % (sys, sys))
            all_aucs=dict()
            max_score, op_scores=get_scores()
            print "max score is ", max_score
            all_active=[]
            all_inactive=[]
            statetracker=[]
            for (i, path) in enumerate(paths.readlines()):
                path=numpy.array(path.split())
                score=numpy.zeros(len(path))
                #h36=numpy.loadtxt('%s/%s%s.h36.dat' % (dir, type, i))
                #in_npxxy=numpy.loadtxt('%s/inactive-npxxy-%s%s.dat' % (dir, type, i))
                #ac_npxxy=numpy.loadtxt('%s/active-npxxy-%s%s.dat' % (dir, type, i))
                #ac_conn=numpy.loadtxt('%s/active-conn-%s%s.dat' % (dir, type, i))
                #in_conn=numpy.loadtxt('%s/inactive-conn-%s%s.dat' % (dir, type, i))
                bulge=numpy.loadtxt('%s/active-h5buldge-%s%s.dat' % (dir, type, i))
                #score=op_path(h36, 'h36', op_scores, score)
                #score=op_path(ac_conn, 'conn', op_scores, score)
                #score=op_path(in_npxxy, 'npxxy', op_scores, score)
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
        #newscores[aucname]=modify_scores(allscores[aucname], newscores[aucname])
    names=['inactive', 'inter', 'active']
    #bar_progress(newscores)
    allprogress=dict()
    limit=False
    for aucname in aucnames:
        #for s in systems:
        #    print "Chi2 Test for %s AUCs and %s Path States" % (aucname, s)
        #    chi2, chi2_s, pval, pval_s, k_chi2, k_pval=get_fried(newscores[aucname][s], names)
        ohandle=open('%s_reduced_allvalues.dat' % aucname, 'w')
        allprogress=dict()
        if aucname=='types':
            limit=True
        else:
            limit=False
        allprogress, modscores=get_aggregate(newscores[aucname], allprogress, names, limit)
        #print "Chi2 Test for %s AUCs and All Path Reduced States" % (aucname)
        #for key in sorted(allprogress.keys()):
        #    for value in allprogress[key]:
        #        ohandle.write('%s\t%s\n' % (key, value)) 
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


