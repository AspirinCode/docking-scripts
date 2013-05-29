import optparse
import pylab
import numpy
import os
from numpy import linalg
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)


def main(type):
    systems=['bi', 'car', 'apo']
    for sys in systems:
        if sys=='bi':
            molecule='Agonist-Bound'
        if sys=='car':
            molecule='Inverse Agonist-Bound'
        if sys=='apo':
            molecule='Apo'
        dir='./%s/structural-data/' % sys
        states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, type), usecols=(0,), dtype=int)
        aucs=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, type), usecols=(1,))
        pops=numpy.loadtxt('%s/Populations.dat' % dir)
        pops=pops[states]
        if sys=='bi':
            cutoff=min(pops)
        else:
            frames=numpy.where(pops>= cutoff)[0]
            pops=pops[frames]
            states=states[frames]
        #names=['path', 'new-path']
        name='path'
        paths=open('%s/%s_paths.txt' % (sys, sys))
        ops=dict()
        orders=['active-h5buldge', 'inactive-npxxy', 'active-npxxy',
                'inactive-conn', 'active-conn', 'h36']
        lims=dict()
        lims['active-h5buldge']=(0,5.0)
        lims['inactive-npxxy']=(0,2.0)
        lims['active-npxxy']=(0,2.0)
        lims['active-conn']=(0,3)
        lims['inactive-conn']=(0,3)
        lims['h36']=(7, 18)
        for (i, path) in enumerate(paths.readlines()):
            path=numpy.array(path.split())
            for order in orders:
                if order not in ops.keys():
                    if order=='h36':
                        ops[order]=numpy.loadtxt('%s/%s%s.h36.dat' % (dir, name, i))
                    else:
                        ops[order]=numpy.loadtxt('%s/%s-%s%s.dat' % (dir, order, name,i))
                else:
                    if order=='h36':
                        ops[order]=numpy.hstack((ops[order],
                            numpy.loadtxt('%s/%s%s.h36.dat' % (dir, name, i))))
                    else:
                        ops[order]=numpy.hstack((ops[order],numpy.loadtxt('%s/%s-%s%s.dat'
                            % (dir, order, name,i))))
            for x in path:
                location=numpy.where(states==int(x))[0]
                if 'aucs' not in ops.keys():
                    ops['aucs']=[]
        dir='./%s/structural-data/' % sys
        #names=['path', 'new-path']
        name='path'
        paths=open('%s/%s_paths.txt' % (sys, sys))
        states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, type), usecols=(0,), dtype=int)
        aucs=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, type), usecols=(1,))
        map=numpy.loadtxt('%s/Mapping.dat' % dir)
        pop=numpy.loadtxt('%s/Populations.dat' % dir)
        ops=dict()
        orders=['active-h5buldge', 'inactive-npxxy', 'active-npxxy',
                'inactive-conn', 'active-conn', 'h36']
        for (i, path) in enumerate(paths.readlines()):
            path=numpy.array(path.split())
            tmp=dict()
            for order in orders:
                if order=='h36':
                    tmp[order]=numpy.loadtxt('%s/%s%s.h36.dat' % (dir, name, i))
                else:
                    tmp[order]=numpy.loadtxt('%s/%s-%s%s.dat' % (dir, order, name,i))
            for (m, x) in enumerate(path):
                location=numpy.where(states==int(x))[0]
                if location.size:
                    if 'aucs' not in ops.keys():
                        ops['aucs']=[]
                        ops['aucs'].append(aucs[location])
                    else:
                        ops['aucs'].append(aucs[location])
                    for order in orders:
                        if order not in ops.keys():
                            ops[order]=[]
                            ops[order].append(tmp[order][m])
                        else:
                            ops[order].append(tmp[order][m])
        print ops.keys()
        for key in ops.keys():
            ops[key]=numpy.array(ops[key])
        ops['aucs']=numpy.array(ops['aucs']).reshape(ops['h36'].shape)
        for x in ops.keys():
            if x is not 'aucs':
                slope, intercept, R, pval, std_err = stats.linregress(ops['aucs'],  ops[x])
                if R**2 > 0.05:
                    if pval< 0.0001:
                        pval=0.0001
                    else:
                        pval=round(pval,4)
                    print sys, x, R, pval
                    pylab.figure()
                    pylab.scatter(ops['aucs'], ops[x], color='k') # , label='%s
                    (ar,br)=polyfit(ops['aucs'], ops[x], 1)
                    xr=polyval([ar,br], ops['aucs'])
                    pylab.plot(ops['aucs'],xr,'r-')
                    pylab.ylim(lims[x])
                    if type=='types':
                        pylab.xlim(0.3, 0.8)
                    else:
                        pylab.xlim(0.5, 1.0)
                    pylab.xlabel('%s %s aucs' % (sys, type))
                    pylab.ylabel(x)
                    pylab.title('%s Simulations R$^2$=%s p=%s' % (molecule,
                        round(R**2,2),pval))
                    pylab.savefig('%s_op%s_%saucs.png' % (sys, x, type) , dpi=300)
                    #pylab.show()


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-t', '--type', dest='type',
                      help='type aucs: agonist, antagonist, types')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(type=options.type)

