import numpy, pylab
systems=['apo', 'bi','car']
type='agonist'
data=dict()
for sys in systems:
    data[sys]=numpy.loadtxt('./%s/%s_types_auc_ci.txt' % (sys, sys), usecols=(1,))

pylab.hist(data['bi'], alpha=0.6, color='r', bins=30, range=(0.2,0.9))
pylab.hist(data['car'], alpha=0.6, color='b', bins=30, range=(0.2,0.9))
pylab.hist(data['apo'], alpha=0.5, color='k', bins=30, range=(0.2,0.9))
pylab.title('discrimination aucs')
pylab.show()
