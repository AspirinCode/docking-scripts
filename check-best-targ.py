#!/bin/python

import numpy, pylab, glob
import optparse
import operator

def main(sys, type ):
    ligand=type
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    best_aucs=dict()
    ci_best_aucs=dict()
    print "on %s" % ligand
    x_aucs=dict()
    x_hici=dict()
    x_lowci=dict()
    xtal=['2rh1', '3p0g']
    alllowci=dict()
    allhici=dict()
    allaucs=dict()
    for state in xtal:
        mfile='../xtal-docking/%s/new-matlab-%s-%s-aucs-95ci.dat' % (ligand, state, ligand)
        value=numpy.loadtxt(mfile, usecols=(0,))
        x_aucs[state]=value
        data=numpy.loadtxt(mfile, usecols=(1,))
        x_lowci[state]=data
        data=numpy.loadtxt(mfile, usecols=(2,))
        x_hici[state]=data
        print "%s %s AUC %s" % (state, ligand, x_hici[state])
    for (n, path) in enumerate(fhandle.readlines()):
        flux=fluxes[n]/fluxes[0]
        states=path.split()
        aucs=[]
        lowci=[]
        hici=[]
        for state in states:
            state=int(state)
            mfile='./%s/targ-matlab-%s-%s-%s-aucs-95ci.dat' % (ligand, sys, state, ligand)
            value=numpy.loadtxt(mfile, usecols=(0,))
            aucs.append(value)
            allaucs[state]=value
            data=numpy.loadtxt(mfile, usecols=(1,))
            lowci.append(data)
            alllowci[state]=data
            data=numpy.loadtxt(mfile, usecols=(2,))
            hici.append(data)
            allhici[state]=data
        # have states in aucs, xtal in x_aucs
        for x in allaucs.keys():
            if type=='agonist':
                state='3p0g'
            if type=='antagonist':
                state='2rh1'
            if allaucs[x] > x_aucs[state]:
                if x in best_aucs.keys():
                    if allaucs[x]>best_aucs[x]:
                        best_aucs[x]=allaucs[x]
                else:
                    best_aucs[x]=allaucs[x]
            if alllowci[x] > x_hici[state]:
                if x in ci_best_aucs.keys():
                    if allaucs[x]>ci_best_aucs[x]:
                        ci_best_aucs[x]=alllowci[x]
                else:
                    ci_best_aucs[x]=alllowci[x]
    ohandle=open('%s_targ_auc_ci.txt' % sys, 'w')
    for key in sorted(allaucs.keys()):
        ohandle.write('%s\t%s\t%s\t%s\n' % (key, allaucs[key], alllowci[key],
            allhici[key]))
    print "CI best %s AUCs" % ligand
    numpy.savetxt('%s_targ_ci_best_aucs.txt' % sys, [i for i in set(ci_best_aucs.keys())],
            fmt='%i')
    print len(ci_best_aucs.keys()), set(ci_best_aucs.keys()), [ci_best_aucs[j]
            for j in ci_best_aucs.keys()]
    print "best %s AUCs" % ligand
    print len(best_aucs.keys()), set(best_aucs.keys()), [best_aucs[j] for j in
            best_aucs.keys()]
    numpy.savetxt('%s_targ_best_aucs.txt' % sys, [i for i in set(best_aucs.keys())], fmt='%i')
    
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    parser.add_option('-t', '--type', dest='type',
            help='type: agonist, antagonist')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, type=options.type)
 
