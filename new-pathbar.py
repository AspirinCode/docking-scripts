#!/bin/python

import numpy, pylab, glob
import optparse

def main(sys, type):
    if type=='agonist':
        other='antagonist'
    elif type=='antagonist':
        other='agonist'
    colors=['magenta', 'cyan']
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    for (n, path) in enumerate(fhandle.readlines()):
        aucs=[]
        lowci=[]
        hici=[]
        pylab.figure()
        print "path %s" % path
        flux=fluxes[n]/fluxes[0]
        states=path.split()
        for state in states:
            state=int(state)
            mfile='./%s/new-matlab-%s-%s-%s-aucs-95ci.dat' % (type, sys, state, type)
            value=numpy.loadtxt(mfile, usecols=(0,))
            aucs.append(value)
            data=numpy.loadtxt(mfile, usecols=(1,))
            lowci.append(value-data)
            data=numpy.loadtxt(mfile, usecols=(2,))
            hici.append(data-value)
        inds=range(0, len(path.split()))
        width=0.8
        pylab.bar(left=inds, height=aucs,  alpha=0.5, width=width,
                color=colors[0], yerr=[lowci, hici], ecolor='k',
                label=type)
        pylab.hold(True)
        pylab.ylim(0.5, 1.0)
        pylab.xticks(numpy.arange(len(states)), states)
        pylab.xlabel('MSM State Used for Docking')
        pylab.ylabel('AUC')
        pylab.title('%s AUCs to Relative Flux %0.2f Path %s' % (type, flux, n))
        pylab.legend(bbox_to_anchor=(1.12, 1.02))
        pylab.savefig('./%s/%s-%s-flux%0.2f-path%s-auc.png' % (type, sys,
            type, flux, n), dpi=300)
        pylab.show()
    
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
 
