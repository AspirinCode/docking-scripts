import numpy, pylab
systems=['bi','car', 'apo']
data=dict()
upper=dict()
lower=dict()
location=dict()
avginactive=dict()
avgactive=dict()
avginterm=dict()
for sys in systems:
    aucs=dict()
    location['active']=[]
    location['inactive']=[]
    location['interm']=[]
    active=numpy.loadtxt('./%s/new-active-states.txt' % (sys), dtype=int)
    inactive=numpy.loadtxt('./%s/new-inactive-states.txt' % (sys), dtype=int)
    states=numpy.loadtxt('./%s/%s_types_auc_ci.txt' % (sys, sys), usecols=(0,),
        dtype=int)
    for (n, x) in enumerate(states):
        if x in active:
            location['active'].append(n)
        elif x in inactive:
            location['inactive'].append(n)
        else:
            location['interm'].append(n)
    location['active']=numpy.array(location['active'])
    location['inactive']=numpy.array(location['inactive'])
    location['interm']=numpy.array(location['interm'])
    data[sys]=numpy.loadtxt('./%s/%s_types_auc_ci.txt' % (sys, sys), usecols=(1,))
    if 'active' not in aucs.keys():
        aucs['active']=data[sys][location['active']]
    else:
        aucs['active']=numpy.hstack((aucs['active'],
            data[sys][location['active']]))
    if 'inactive' not in aucs.keys():
        aucs['inactive']=data[sys][location['inactive']]
    else:
        aucs['inactive']=numpy.hstack((aucs['inactive'],
            data[sys][location['inactive']]))
    if 'interm' not in aucs.keys():
        aucs['interm']=data[sys][location['interm']]
    else:
        aucs['interm']=numpy.hstack((aucs['interm'],
            data[sys][location['interm']]))
    #upper[sys]=numpy.loadtxt('./%s/%s_types_auc_ci.txt' % (sys, sys), usecols=(2,))
    #lower[sys]=numpy.loadtxt('./%s/%s_types_auc_ci.txt' % (sys, sys), usecols=(3,))

    avgactive[sys]=numpy.mean(aucs['active'])
    avginactive[sys]=numpy.mean(aucs['inactive'])
    avginterm[sys]=numpy.mean(aucs['interm'])
    #pylab.hist(aucs['active'], alpha=0.6, color='r', bins=30, range=(0.2,0.9))
    #pylab.hist(aucs['inactive'], alpha=0.6, color='b', bins=30, range=(0.2,0.9))
    #pylab.hist(aucs['interm'], alpha=0.5, color='k', bins=30, range=(0.2,0.9))
    #pylab.title('%s/discrimination aucs' % sys)

types=['ro', 'rx', 'r-', 'bo', 'bx', 'b-', 'ko', 'kx', 'k-']
colors=['r', 'b', 'k']
for (n, sys) in enumerate(systems):
    pylab.plot([n,], avgactive[sys], 'o', color=colors[n], label='%s active' % sys) 
    pylab.plot([n,], avginactive[sys], 'x', color=colors[n],  label='%s inactive' % sys) 
    pylab.plot([n,], avginterm[sys], '+', color=colors[n], label='%s intermediate' % sys) 
pylab.xticks(range(0, len(systems)), ['agonist', 'inv. agonist', 'apo'])
pylab.xlim(-1,3)
pylab.ylim(0.4, 0.7)
pylab.show()
