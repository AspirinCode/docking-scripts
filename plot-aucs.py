import numpy, pylab
systems=['apo', 'bi','car']
type='agonist'
types=['agonist', 'antagonist']
for type in types:
    data=dict()
    pops=dict()
    for sys in systems:
        data[sys]=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, type), usecols=(1,))
        states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, type),
                usecols=(0,), dtype=int)
        p=numpy.loadtxt('./%s/Populations.dat' % sys)
        pops[sys]=p[states]
    cmap=pylab.get_cmap("RdYlBu")
    pylab.figure()
    pylab.scatter([-numpy.log(i) for i in pops['bi']], data['bi'],  color='r',alpha=0.6)
    pylab.scatter([-numpy.log(i) for i in pops['car']], data['car'], color='b', alpha=0.6)
    pylab.scatter([-numpy.log(i) for i in pops['apo']], data['apo'], color='k', alpha=0.5)
    pylab.xlabel('population')
    pylab.ylabel('aucs')
    pylab.title('%s aucs population' % type)
    pylab.figure()
    pylab.hist(data['bi'], color='r', alpha=0.6, bins=20, range=(0.5,1.0))
    pylab.hist(data['car'], color='b', alpha=0.6, bins=20, range=(0.5,1.0))
    pylab.hist(data['apo'], color='k', alpha=0.5, bins=20, range=(0.5,1.0))
    pylab.title('%s aucs' % type)
    pylab.show()
