import numpy, pylab, random
from scipy import stats
ref=numpy.loadtxt('bi_agonist_aucs.txt', usecols=(1,))
rand=numpy.loadtxt('bi_rand_aucs.txt', usecols=(1,))
frames=range(0, len(ref))
random.shuffle(frames)
frames=frames[:len(rand)]
ref=ref[frames]
print stats.chisquare(ref, rand)
pylab.hist(ref, bins=15, range=[0.7,0.95], alpha=0.6, label='MSM States', normed=True)
pylab.hold(True)
pylab.hist(rand, bins=15, range=[0.7,0.95], alpha=0.6, label='Random States', normed=True)
pylab.legend()
pylab.xlabel('Docking Performance AUC')
pylab.ylabel('Normed Probability')
pylab.show()
pylab.savefig('compare-rand.png', dpi=300)
