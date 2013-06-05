import optparse
import numpy, pylab, random
from math import  *
from contingency import chi2_contingency
import scipy.stats as stat

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)


def ztest(teststat, nullvalue,  se, alpha =0.05, side=0):
    """
    Normal test of a sample statistic.
    Return:
       pvalue
       Z test statistic
       critical value(s)

    Arguments:
      teststat- test statistic
      se -    standard error of sample statistic
      alpha - significance level
      side  - -1 left-sided test
             0  double sided test
             +1 right-sided test
    """
    pnorm = stat.norm.cdf
    qnorm = stat.norm.ppf
    Ztest = (teststat-nullvalue)/se
    print "Ztest=",  Ztest
    if side ==0:
       pvalue = pnorm(Ztest)
       if Ztest > 0.0:
           pvalue = 1.0 -pvalue
       pvalue *= 2.0
       zcrit1 = qnorm(alpha/2)
       zcrit2 = qnorm(1-alpha/2.0)
       return pvalue, Ztest, (zcrit1,zcrit2)
    elif side == -1:
       pvalue = pnorm(Ztest)
       zcrit = qnorm(alpha)
       return pvalue, Ztest, zcrit
    else:
       pvalue = 1- pnorm(Ztest)
       zcrit  = qnorm(1.0-alpha)
       return pvalue, Ztest, zcrit

def twosampleproptest(x1,  n1,   x2,  n2,  alpha,  side):
    p1hat = float(x1)/n1
    p2hat = float(x2)/n2
    phat = float(x1 +x2)/(n1+n2)
    print "p1hat, p2hat",  p1hat,  p2hat
    samplestat =  p1hat - p2hat
    se = sqrt(phat*(1.0-phat) * (1.0/n1 + 1.0/n2))
    print "se %s" % se
    print "%s-sided" % side
    return ztest(p2hat, p1hat,  se,  alpha,  side)

def check(data, xtal):
    count=0
    for x in data:
        if x > xtal:
            count+=1
    return count

def main(sys, type):
    xtal_low=numpy.loadtxt('../xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' %
        (type, type), usecols=(1,))
    xtal_hi=numpy.loadtxt('../xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' %
        (type, type) ,usecols=(2,))
    xtal_auc=numpy.loadtxt('../xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' %
        (type, type),
            usecols=(0,))
    ref=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(1,))
    ref_low=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(2,))
    ref_hi=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(3,))
    rand=numpy.loadtxt('%s_rando_%s_auc_ci.txt' % (sys, type), usecols=(1,))
    rand=rand[:100]
    rand_low=numpy.loadtxt('%s_rando_%s_auc_ci.txt' % (sys, type), usecols=(2,))
    rand_low=rand_low[:100]
    rand_hi=numpy.loadtxt('%s_rando_%s_auc_ci.txt' % (sys, type), usecols=(3,))
    rand_hi=rand_hi[:100]
    t, pval=stat.stats.ttest_ind(ref, rand)
    print "t= ", t, "pval= ", pval
    chi2, pval=stat.stats.chisquare(ref, rand)
    print "chi2= ", chi2, "pval= ", pval
    count=check(rand, xtal_auc)
    count2=check(rand_low, xtal_hi)
    print "%s Random State AUCs > Xtal AUC" % count
    count=check(ref, xtal_auc)
    count2=check(ref_low, xtal_hi)
    print "%s MSM State AUCs > Xtal AUC" % count
    pylab.figure()
    pylab.hold(True)
    test2=pylab.hist(rand, bins=15, range=[0.7,0.95], color='cyan', label='Random States', normed=True)
    test=pylab.hist(ref, bins=15, range=[0.7,0.95], alpha=0.6, color='magenta',  label='MSM States', normed=True)
    maxval=numpy.hstack((test[0], test2[0])).max()
    pylab.plot([xtal_hi]*(int(maxval)+2), range(0, int(maxval)+2),
            linestyle='dotted', color='k', label='xtal upper CI') 
    pylab.ylim(0, maxval+1)
    pylab.legend()
    pylab.xlabel('Docking Performance AUC')
    pylab.ylabel('Normed Probability')
    pylab.savefig('compare-rand.png', dpi=300)
    #pylab.show()
    
    

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, type=options.type)

    
