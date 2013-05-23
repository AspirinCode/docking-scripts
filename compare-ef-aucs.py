import numpy, pylab
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

auc_state=numpy.loadtxt('bi_agonist_aucs.txt', usecols=(0,), dtype=int)
aucs=numpy.loadtxt('bi_agonist_aucs.txt', usecols=(1,))
order=numpy.argsort(aucs)
sort_auc_state=auc_state[order][::-1]
sort_aucs=aucs[order][::-1]

# ref is EF
ef_state=numpy.loadtxt('allstate_10ef.dat', usecols=(0,), dtype=int)
efs=numpy.loadtxt('allstate_10ef.dat', usecols=(1,))
order=numpy.argsort(efs)
sort_ef_state=ef_state[order][::-1]
sort_efs=efs[order][::-1]


new_aucs=[]
new_auc_states=[]
new_efs=[]
new_ef_states=[]
for state in sort_ef_state:
    index=numpy.where(sort_auc_state==state)[0]
    if index.size:
        new_auc_states.append(state)
        new_aucs.append(sort_aucs[index])
        new_ef_states.append(state)
        index=numpy.where(sort_ef_state==state)[0]
        new_efs.append(sort_efs[index])

new_aucs=[x[0] for x in new_aucs]
new_efs=[x[0] for x in new_efs]
print "EF", new_ef_states[:10]
print "AUC", new_auc_states[:10]
slope, intercept, R, pval, std_err = stats.linregress(new_aucs, new_efs)
print R, pval
pylab.scatter(new_aucs, new_efs)
(ar,br)=polyfit(new_aucs, new_efs, 1)
xr=polyval([ar,br],new_aucs)
pylab.plot(new_aucs,xr,'r-')
pylab.xlabel('AUCs')
pylab.ylabel('EF at 10%% DB')
pylab.title('R$^2$=%s, pval=%s' % (round(R**2,2), round(pval, 2)))
pylab.savefig('compare-ef10-aucs.png', dpi=300)
