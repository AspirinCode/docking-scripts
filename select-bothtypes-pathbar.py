#!/bin/python
import numpy, pylab, glob
import optparse

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)



def main(sys, type, state):
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
        if state in path.split():
            print "path %s" % n
            flux=fluxes[n]/fluxes[0]
            states=path.split()
            pylab.figure()
            for (m, ligand) in enumerate(ligands):
                print "on %s" % ligand
                aucs=[]
                lowci=[]
                hici=[]
                for state in states:
                    state=int(state)
                    mfile='./%s/new-matlab-%s-%s-%s-aucs-95ci.dat' % (ligand, sys, state, ligand)
                    value=numpy.loadtxt(mfile, usecols=(0,))
                    aucs.append(value)
                    data=numpy.loadtxt(mfile, usecols=(1,))
                    lowci.append(value-data)
                    data=numpy.loadtxt(mfile, usecols=(2,))
                    hici.append(data-value)
                x_aucs=dict()
                x_hici=dict()
                x_lowci=dict()
                xtal=['2rh1', '3p0g']
                for state in xtal:
                    mfile='../xtal-docking/%s/new-matlab-%s-%s-aucs-95ci.dat' % (ligand, state, ligand)
                    value=numpy.loadtxt(mfile, usecols=(0,))
                    x_aucs[state]=value
                    data=numpy.loadtxt(mfile, usecols=(1,))
                    x_lowci[state]=value-data
                    data=numpy.loadtxt(mfile, usecols=(2,))
                    x_hici[state]= data-value
                total_aucs=numpy.hstack((numpy.array(x_aucs['3p0g']), aucs,numpy.array(x_aucs['2rh1'])))
                total_lowci=[float(x_lowci['3p0g'])]+lowci+[float(x_lowci['2rh1'])]
                total_hici=[float(x_hici['3p0g'])]+hici+[float(x_hici['2rh1'])]
                print 'agonist', total_lowci, total_hici
                newstates=['3p0g'] + [str(i) for i in states] + ['2rh1']
                inds=range(0, len(newstates))
                width=0.8
                pylab.bar(left=inds, height=total_aucs,  alpha=0.5, width=width,
                        color=colors[m], yerr=[total_lowci, total_hici], ecolor='k',
                        label=ligand)
                pylab.hold(True)
            pylab.ylim(0.5, 1.0)
            pylab.xticks(numpy.arange(len(newstates)), newstates)
            pylab.xlabel('MSM State Used for Docking')
            pylab.ylabel('AUC')
            if sys=='bi':
                name='Agonist-bound'
            elif sys=='car':
                name='Antagonist-bound'
            elif sys=='apo':
                name='Apo'
            pylab.title('AUCs to %s Path %s %s%% Flux' % (name, n, int(flux*100)))
            pylab.legend(bbox_to_anchor=(1.12, 1.02))
            pylab.savefig('./select-%s-flux%0.2f-path%s-auc.png' % (sys, flux, n), dpi=300)
            pylab.show()
            break
    
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    parser.add_option('-t', '--type', dest='type',
            help='type: agonist, antagonist')
    parser.add_option('-k', '--state', dest='state',
            help='state to find path to plot')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, type=options.type, state=options.state)

