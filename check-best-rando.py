#!/bin/python

import numpy, pylab, glob
import optparse
import operator

def main(sys, type ):
    ligand=type
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    aucs=dict()
    best_aucs=dict()
    ci_best_aucs=dict()
    print "on %s" % ligand
    x_aucs=dict()
    x_hici=dict()
    x_lowci=dict()
    xtal=['2rh1', '3p0g']
    for state in xtal:
        mfile='../xtal-docking/%s/new-matlab-%s-%s-aucs-95ci.dat' % (ligand, state, ligand)
        value=numpy.loadtxt(mfile, usecols=(0,))
        x_aucs[state]=value
        data=numpy.loadtxt(mfile, usecols=(1,))
        x_lowci[state]=data
        data=numpy.loadtxt(mfile, usecols=(2,))
        x_hici[state]=data
        print "%s %s AUC %s" % (state, ligand, x_hici[state])
    aucs=dict()
    lowci=dict()
    hici=dict()
    for state in range(1,10):
        state=int(state)
        mfile='./%s/new-matlab-rand-%s-%s-%s-aucs-95ci.dat' % (ligand, sys, state, ligand)
        value=numpy.loadtxt(mfile, usecols=(0,))
        aucs[state]=value
        data=numpy.loadtxt(mfile, usecols=(1,))
        lowci[state]=data
        data=numpy.loadtxt(mfile, usecols=(2,))
        hici[state]=data
    for x in sorted(aucs.keys()):
        if type=='agonist':
            state='3p0g'
        if type=='antagonist':
            state='2rh1'
        if aucs[x] > x_aucs[state]:
            best_aucs[x]=aucs[x]
        else:
            pass
        if lowci[x] > x_hici[state]:
            ci_best_aucs[x]=aucs[x]
        else:
            pass
    print "CI best %s AUCs" % len(ci_best_aucs.keys())
    print "not-CI best %s AUCs" % len(best_aucs.keys())
    #print sorted(ci_best_aucs.iteritems(),key=operator.itemgetter(1))[::-1]
    ohandle=open('%s_rand_aucs.txt' % sys, 'w')
    for x in sorted(aucs.keys()):
        ohandle.write('%s\t%s\n' % (x, aucs[x]))
    #best=[]
    #for j in set([x[1][1] for x in sorted(ci_best_aucs.iteritems(),
    #    key=operator.itemgetter(1))]):
    #    best.append(int(j))
    #numpy.savetxt('ci_best_aucs.txt', numpy.array(best), fmt='%i')
    
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
 
