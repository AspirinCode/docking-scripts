#!/bin/python

import numpy, pylab, glob
import optparse
import operator

def main(sys ):
    ligand='eval-types'
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    high_aucs=dict()
    low_aucs=dict()
    ci_high_aucs=dict()
    ci_low_aucs=dict()
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
            mfile='./%s/new-matlab-%s-%s-%s-aucs-95ci.dat' % (ligand, sys, state, ligand)
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
            if allaucs[x] > x_aucs['3p0g']:
                if x in high_aucs.keys():
                    if allaucs[x]>high_aucs[x]:
                        high_aucs[x]=allaucs[x]
                else:
                    high_aucs[x]=allaucs[x]
            if alllowci[x] > x_hici['3p0g']:
                if x in ci_high_aucs.keys():
                    if allaucs[x]>ci_high_aucs[x]:
                        ci_high_aucs[x]=alllowci[x]
                else:
                    ci_high_aucs[x]=alllowci[x]

            if allaucs[x] < x_aucs['2rh1']:
                if x in low_aucs.keys():
                    if allaucs[x]<low_aucs[x]:
                        low_aucs[x]=allaucs[x]
                else:
                    low_aucs[x]=allaucs[x]
            if allhici[x] < x_lowci['2rh1']:
                if x in ci_low_aucs.keys():
                    if allaucs[x]<ci_low_aucs[x]:
                        ci_low_aucs[x]=allhici[x]
                else:
                    ci_low_aucs[x]=allhici[x]

    ohandle=open('%s_types_auc_ci.txt' % sys, 'w')
    for key in sorted(allaucs.keys()):
        ohandle.write('%s\t%s\t%s\t%s\n' % (key, allaucs[key], alllowci[key],
            allhici[key]))
    print "CI high %s AUCs" % ligand
    numpy.savetxt('%s_types_ci_high_aucs.txt' % sys, [i for i in set(ci_high_aucs.keys())],
            fmt='%i')
    print len(ci_high_aucs.keys()), set(ci_high_aucs.keys()), [ci_high_aucs[j]
            for j in ci_high_aucs.keys()]
    print "CI low %s AUCs" % ligand
    numpy.savetxt('%s_types_ci_low_aucs.txt' % sys, [i for i in set(ci_low_aucs.keys())],
            fmt='%i')
    print len(ci_low_aucs.keys()), set(ci_low_aucs.keys()), [ci_low_aucs[j]
            for j in ci_low_aucs.keys()]
    print "high %s AUCs" % ligand
    print len(high_aucs.keys()), set(high_aucs.keys()), [high_aucs[j] for j in
            high_aucs.keys()]
    numpy.savetxt('%s_types_high_aucs.txt' % sys, [i for i in set(high_aucs.keys())], fmt='%i')
    print "low %s AUCs" % ligand
    print len(low_aucs.keys()), set(low_aucs.keys()), [low_aucs[j] for j in low_aucs.keys()]
    numpy.savetxt('%s_types_low_aucs.txt' % sys, [i for i in set(low_aucs.keys())], fmt='%i')
    
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys)
 
