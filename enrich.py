#!/bin/python
import os
import numpy, pylab, glob
import operator 
import random
import optparse

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)

def get_best(d):
    sorted_data=[ value for value in sorted(d.iteritems(),
        key=operator.itemgetter(1))][::-1]
    vals=numpy.array([x[1] for x in sorted_data])
    keys=numpy.array([x[0] for x in sorted_data])
    best=[]
    for (c, i) in enumerate(vals):
        if (i-max(vals))/max(vals) < 0.1:
            best.append(sorted_data[c])
        else:
            break
    return best

def get_data(type, sys, state):
    data=dict()
    ligands='%s/ranked-%s-aln-state%s.cent-ligand%s.out' % ( type, sys, state, type)
    decoys='%s/ranked-%s-aln-state%s.cent-ligand%s-decoys.out' % ( type, sys, state, type)
    data_l=numpy.loadtxt(ligands, usecols=(1,))
    if not os.path.exists(decoys):
        return 0, 0, 0, 0, 0, 0, 0
    else:
        data_d=numpy.loadtxt(decoys, usecols=(1,))
        total=len(data_d)+len(data_l)
        total_hits=len(data_l)
        for (n,x) in enumerate(data_l):
            name='hit%s' % n
            data[name]=x
        for (n,x) in enumerate(data_d):
            name='decoy%s' % n
            data[name]=x
        # ef(subset size)=(hit_s/n_s)/(hit_t/n_t)
        sorted_data=numpy.array(sorted(data.iteritems(),
            key=operator.itemgetter(1)))
        hitcount=0
        true=float(len(data_l))/(len(data_l) + len(data_d))*100
        ef=[]
        ligands=[]
        for (n, x) in enumerate(sorted_data[::-1]):
            if 'hit' in x[0]:
                hitcount+=1
            ef.append((float(hitcount)/(n+1))/(float(total_hits)/total))
            ligands.append(100*(float(hitcount)/total_hits))
        percent=[100*(float(j)/total) for j in numpy.arange(0, total)]
        ef=numpy.array(ef)
        return total_hits, data, numpy.array(ligands), numpy.array(ef), numpy.array(percent), true, total

def plot_ef(sys, total_hits, data, percent, ef, ligands, state, total):
    pylab.scatter([numpy.log10(j) for j in percent], ligands, color='red',
            label='state%s' % (state))
    #pylab.scatter(percent, ef, color='red')
    hitcount=0
    ligands=[]
    for (n, x) in enumerate(data.keys()):
        if 'hit' in x:
            hitcount+=1
        ligands.append(100*(float(hitcount)/total_hits))
    pylab.scatter([numpy.log10(j) for j in percent], ligands, color='k',
            label='random')
    pylab.hold(True)
    cutoff=0.1
    num=int(total*cutoff)
    frames=numpy.arange(0, num)
    pylab.plot([cutoff*10]*100, range(0,100), '--', color='k', label='EF %s%%' % (cutoff*100))
    pylab.title('Enrichment State %s' % state)
    pylab.xlim(-1, 2) # x log goes from -1 (0.1) to 2 (100)
    pylab.xticks(numpy.arange(-2.0, 3.0, 1.0), ['0', '0.1', '1.0', '10', '100'])
    pylab.ylim(0,100)
    pylab.legend(loc=2)
    pylab.xlabel('% of Database Screened')
    pylab.ylabel('% of Hits Found')
    pylab.savefig('%s_state%s_enrichment.png' % (sys, state) , dpi=300)
    #pylab.show()

def main(sys, type):
    file='%s_states.txt' % sys
    states=numpy.loadtxt(file)
    ef_2=dict()
    ef_5=dict()
    ef_10=dict()
    ef_max=dict()
    for state in states:
        state=int(state)
        print "on state %s" % state
        total_hits, data, ligands, ef, percent, true, total=get_data(type, sys, state)
        if total_hits==0:
            continue
        else:
            print "random enrich is %s" % (total_hits/float(total))
            index=numpy.where(ef==max(ef))[0]
            if len(index)>1:
                frames=numpy.arange(0, index[-1])
                fraction=float(frames[-1])/len(ef)
            else:
                fraction=float(index)/len(ef)
            ef_max['percent']=round(fraction,2)
            print "max EF %s at" % max(ef), round(fraction,2)
            ef_max['val']=max(ef)
            for (d, cutoff) in zip([ef_2, ef_5, ef_10], [0.02, 0.05, 0.1]):
                num=int(total*cutoff)
                frames=numpy.arange(0, num)
                d[state]=ef[frames][-1]
            # plot percent ligands recovered with log percent
            #plot_ef(sys, total_hits, data, percent, ef, ligands, state, total)
    count=0
    colors=['magenta', 'cyan', 'grey']
    for (d, cutoff) in zip([ef_2, ef_5, ef_10], [0.02, 0.05, 0.1]):
        tmp=[]
        best=get_best(d)
        print cutoff, best
        ohandle=open('allstate_%sef.dat' % int(cutoff*100), 'w')
        for key in sorted(d.keys()):
            ohandle.write('%s\t%s\n' % (key, d[key]))
            tmp.append(d[key])
        pylab.hist(tmp, label='%s%% DB Screened' % int(cutoff*100), alpha=0.6,
                color=colors[count], bins=24, range=[0,7], normed=True)
        pylab.hold(True)
        pair=best[0]
        #pylab.plot([pair[1]]*100, numpy.arange(0,1, 0.01), color=colors[count],
        #        label='%s' % pair[0])
        count+=1
    pylab.xlabel('EF Values')
    pylab.legend()
    pylab.ylabel('Normalized Probability')
    pylab.savefig('%s_ef_values.png' % sys, dpi=300)
            
    
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--sys', dest='sys',
                      help='system')
    parser.add_option('-t', '--type', dest='type',
                      help='type ligands')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, type=options.type)

