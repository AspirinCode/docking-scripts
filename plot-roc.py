#!/bin/python

import numpy, pylab, glob
import optparse

def main(sys, type ):
    if type=='agonist':
        other='antagonist'
    elif type=='antagonist':
        other='agonist'
    ligands=[other, type]
    colors=['blue', 'red']
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    for (n, path) in enumerate(fhandle.readlines()):
        print "path %s" % n
        flux=fluxes[n]/fluxes[0]
        states=path.split()
        pylab.figure()
        for (m, ligand) in enumerate(ligands):
            print "on %s" % ligand
            count=0
            for state in states:
                state=int(state)
                mfile='./%s/new-matlab-%s-%s-%s-roc.dat' % (ligand, sys, state, ligand)
                x=numpy.loadtxt(mfile, usecols=(1,))
                y=numpy.loadtxt(mfile, usecols=(0,))
                if count==0:
                    pylab.plot(x,y, label=ligand,
                        color=colors[m])
                else:
                    pylab.plot(x,y, color=colors[m])
                count+=1
            xtal=['3p0g', '2rh1']
            xcolors=['#000066', '#660033']
            for state in xtal:
                mfile='../xtal-docking/%s/new-matlab-%s-%s-roc.dat' % (ligand, state, ligand)
                x=numpy.loadtxt(mfile, usecols=(1,))
                y=numpy.loadtxt(mfile, usecols=(0,))
                if state=='3p0g':
                    pylab.plot(x,y, '--', label='%s\t%s' % (state, ligand),
                        color=xcolors[m], linewidth=3)
                else:
                    pylab.plot(x,y, '-', label='%s\t%s' % (state, ligand),
                        color=xcolors[m], linewidth=3)
        pylab.xlabel('false positive rate, 1-specificity')
        pylab.ylabel('true positive rate, sensitivity')
        pylab.ylim(0,1)
        if sys=='bi':
            name='Agonist-bound'
        elif sys=='car':
            name='Antagonist-bound'
        elif sys=='apo':
            name='Apo'
        pylab.title('ROC for %s Receptor Path %s Relative Flux %0.2f' % (name, n, flux))
        pylab.legend(bbox_to_anchor=(1.12, 0.8))
        pylab.savefig('./top-%s-flux%0.2f-path%s-roc.png' % (sys, flux, n), dpi=300)
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
 
