import optparse
from DockingStats import *
import numpy, pylab, random
from math import  *
from contingency import chi2_contingency
import scipy.stats as stat

font = {'family' : 'sans-serif',
                'sans-serif':['Helvetica'],
                                        'size'   : 14}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 10,
                  'legend.linewidth': 2}
pylab.rcParams.update(params)



def main(sys, type, refname, testname):
    testlabel=get_label(testname)
    reflabel=get_label(refname)
    states=numpy.loadtxt('%s/%s_%s_%s_auc_ci.txt' % (sys, sys, refname, type), usecols=(0,))
    filter=numpy.loadtxt('%s/filter_states.txt' % sys)
    new=[]
    for (n, x) in enumerate(states):
        if float(x) in filter:
            new.append(n)
    new=numpy.array(new)

    data=dict()
    data['reference']=dict()
    data['test']=dict()
    data['xtal']=dict()
    names=['%s/%s_%s_%s_auc_ci.txt' % (sys, sys, refname, type),
            '%s/%s_%s_%s_auc_ci.txt' % (sys, sys, testname, type), 'xtal-docking/%s/new-matlab-3p0g-%s-aucs-95ci.dat' % (type, type)]  
    for (key, name) in zip(['reference', 'test', 'xtal'], names):
        d=data[key]
        if 'xtal' in name:
            n=0
            filter=False
        else:
            filter=True
            n=1
        d['avg']=numpy.loadtxt(name, usecols=(n,))
        if filter==True:
            d['avg']=d['avg'][new]
        n+=1
        d['low']=numpy.loadtxt(name, usecols=(n,))
        if filter==True:
            d['low']=d['low'][new]
        n+=1
        d['high']=numpy.loadtxt(name, usecols=(n,))
        if filter==True:
            d['high']=d['high'][new]

    print "%s mean= " % testlabel, numpy.mean(data['test']['avg']), "%s mean= " % reflabel, numpy.mean(data['reference']['avg'])
    print "%s variance= " % testlabel, numpy.var(data['test']['avg']), "%s variance= " % reflabel, numpy.var(data['reference']['avg'])
    #statistical tests of difference in population
    t, pval=ttest_ind(data['test']['avg'], data['reference']['avg'], axis=0, equal_var=False)
    print "population difference t= ", t, "pval= ", pval
    chi2, pval=stat.stats.chisquare(data['reference']['avg'], data['test']['avg'])
    print "population difference chi2= ", chi2, "pval= ", pval
    
    #statistical tests of difference in proportion states < xtal
    counts=dict()
    counts['test']=dict()
    counts['reference']=dict()
    for key in counts.keys():
        for level in data['test'].keys():
            if level=='high':
                #counts[key][level]=check(data[key]['avg'], data['xtal'][level]) #constant xtal
                counts[key][level]=check(data[key][level], data['xtal']['low']) #opposite level for xtal
                #counts[key][level]=check(data[key][level], data['xtal']['avg']) # constant avg xtal
                print "Freq. AUCs > Xtal %s loose criteria: " % key, counts[key][level]
            elif level=='low':
                #counts[key][level]=check(data[key]['avg'], data['xtal'][level]) #opposite level for xtal
                counts[key][level]=check(data[key][level], data['xtal']['high']) # strict
                #counts[key][level]=check(data[key][level], data['xtal']['avg']) # strict
                print "Freq. AUCs > Xtal %s strict criteria: " % key, counts[key][level]
            elif level=='avg':
                #counts[key][level]=check(data[key]['avg'], data['xtal'][level]) #opposite level for xtal
                counts[key][level]=check(data[key][level], data['xtal'][level]) 
                #counts[key][level]=check(data[key][level], data['xtal']['avg']) 
                print "Freq. AUCs > Xtal %s average criteria: " % key, counts[key][level]
    side=1
    print "%s-Sided Test" % side
    for (level, c) in zip(['strict', 'avg', 'loose'], ['low', 'avg', 'high']):
        pvalue, Ztest, zcrit=twosampleproptest(counts['test'][c], len(data['test'][c]),
                counts['reference'][c], len(data['reference'][c]),
                alpha=0.05,  side=side)
        print "%s Z= " % level, Ztest, "zcrit= ", zcrit, "pvalue= ", pvalue
    pylab.figure()
    pylab.hold(True)
    if type=='types':
        binrange=[0.2, 1.00]
        binum=20
    else:
        binrange=[0.6, 1.00]
        binum=15
    test2=pylab.hist(data['test']['avg'], bins=binum, range=binrange,
            color='cyan', label='%s States' % testlabel, normed=True)
    test=pylab.hist(data['reference']['avg'], bins=binum, range=binrange, alpha=0.6, color='magenta',  label='MSM States', normed=True)
    maxval=numpy.hstack((test[0], test2[0])).max()
    for x in ['low', 'avg', 'high']:
        if x=='avg':
            format='k-'
            label='Xtal'
            pylab.plot([data['xtal'][x]]*(int(maxval)+2), range(0, int(maxval)+2), format, label=label) 
        else:
            format='k--'
        pylab.plot([data['xtal'][x]]*(int(maxval)+2), range(0, int(maxval)+2), format)
    if type=='types':
        pylab.plot([0.5]*(int(maxval)+2), range(0, int(maxval)+2), 'k-.')
    pylab.ylim(0, maxval)
    pylab.legend()
    pylab.xlabel('%s AUC' % get_label(type))
    pylab.ylabel('Normed Probability')
    pylab.title('%s %s vs. %s States' % (get_label(sys), reflabel, testlabel))
    pylab.savefig('compare-%s-%s-rand.png' % (sys, type), dpi=300)
    pylab.show()
    
    

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    parser.add_option('-r', '--refname', dest='refname',
                      help='reference AUC set')
    parser.add_option('-x', '--testname', dest='testname',
                      help='test AUC set')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, type=options.type, refname=options.refname,
            testname=options.testname)

    
