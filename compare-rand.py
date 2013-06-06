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

def ttest_finish(df,t):
    """Common code between all 3 t-test functions."""
    prob = stat.distributions.t.sf(numpy.abs(t), df) * 2  # use numpy.abs to get upper tail
    if t.ndim == 0:
        t = t[()]

    return t, prob

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
    #print "Ztest=",  Ztest
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

def ttest_ind(a, b, axis=0, equal_var=True):
    a, b, axis = stat._support._chk2_asarray(a, b, axis)
    v1 = numpy.var(a, axis, ddof=1)
    v2 = numpy.var(b, axis, ddof=1)
    n1 = a.shape[axis]
    n2 = b.shape[axis]

    if (equal_var):
        df = n1 + n2 - 2
        svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / float(df)
        denom = numpy.sqrt(svar * (1.0 / n1 + 1.0 / n2))
    else:
        vn1 = v1 / n1
        vn2 = v2 / n2
        df = ((vn1 + vn2)**2) / ((vn1**2) / (n1 - 1) + (vn2**2) / (n2 - 1))

        # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
        # Hence it doesn't matter what df is as long as it's not NaN.
        df = numpy.where(numpy.isnan(df), 1, df)
        denom = numpy.sqrt(vn1 + vn2)

    d = numpy.mean(a, axis) - numpy.mean(b, axis)
    t = numpy.divide(d, denom)
    t, prob = ttest_finish(df, t)

    return t, prob

def twosampleproptest(x1,  n1,   x2,  n2,  alpha,  side):
    p1hat = float(x1)/n1
    p2hat = float(x2)/n2
    phat = float(x1 +x2)/(n1+n2)
    #print "p1hat, p2hat",  p1hat,  p2hat
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
    states=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(0,))
    filter=numpy.loadtxt('filter_states.txt')
    new=[]
    for (n, x) in enumerate(states):
        if float(x) in filter:
            new.append(n)
    new=numpy.array(new)

    reference=dict()
    random=dict()
    xtal_auc=numpy.loadtxt('../xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' % (type, type), usecols=(0,))
    xtal_low=numpy.loadtxt('../xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' % (type, type), usecols=(1,))
    xtal_hi=numpy.loadtxt('../xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' % (type, type), usecols=(2,))
    ref=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(1,))
    ref=ref[new]
    ref_low=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(2,))
    ref_low=ref_low[new]
    ref_hi=numpy.loadtxt('%s_%s_auc_ci.txt' % (sys, type), usecols=(3,))
    ref_hi=ref_hi[new]
    rand=numpy.loadtxt('%s_rando_%s_auc_ci.txt' % (sys, type), usecols=(1,))
    rand=rand[:len(ref)]
    rand_low=numpy.loadtxt('%s_rando_%s_auc_ci.txt' % (sys, type), usecols=(2,))
    rand_low=rand_low[:len(ref)]
    rand_hi=numpy.loadtxt('%s_rando_%s_auc_ci.txt' % (sys, type), usecols=(3,))
    rand_hi=rand_hi[:len(ref)]
    print "random mean= ", numpy.mean(rand), "ref mean= ", numpy.mean(ref)
    print "random variance= ", numpy.var(rand), "ref variance= ", numpy.var(ref)
    #statistical tests of difference in population
    t, pval=ttest_ind(ref, rand, axis=0, equal_var=False)
    print "population difference t= ", t, "pval= ", pval
    chi2, pval=stat.stats.chisquare(ref, rand)
    print "population difference chi2= ", chi2, "pval= ", pval
    
    #statistical tests of difference in proportion states < xtal
    rand_count=check(rand, xtal_auc)
    rand_count_loose=check(rand_hi, xtal_low)
    print "%s Random State AUCs > Xtal AUC" % rand_count
    ref_count=check(ref, xtal_auc)
    ref_count_loose=check(ref_hi, xtal_low)
    ref_count_=check(ref_hi, xtal_low)
    print "%s MSM State AUCs > Xtal AUC" % ref_count
    pvalue, Ztest, zcrit=twosampleproptest(rand_count,  len(rand),   ref_count,
            len(ref),  alpha=0.05,  side=1)
    print "Frequence > Xtal Z= ", Ztest, "zcrit= ", zcrit, "pvalue= ", pvalue
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
    pylab.show()
    
    

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

    
