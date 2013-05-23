#!/bin/python

import numpy, pylab, glob
import optparse
import operator

def main(sys):
    types=['agonist']
    for type in types:
        print "on %s" % type
        files=glob.glob('%s/new-matlab-%s-state*-%s*95*' % (type, sys, type))
        best=dict()
        bestci=dict()
        for file in files:
            state=file.split('-%s' % type)[0].split('-state')[1]
            value=numpy.loadtxt(file, usecols=(0,))
            data1=numpy.loadtxt(file, usecols=(1,))
            data2=numpy.loadtxt(file, usecols=(2,))
            if state in best.keys():
                if value!=best[state]:
                    print "issue"
                    import pdb
                    pdb.set_trace()
            else:
                best[state]=value
                bestci[state]=[(data1, data2)]
        xtal=dict()
        file='./prod2/new-matlab-2rh1-%s-aucs-95ci.dat' % type
        value=numpy.loadtxt(file, usecols=(0,))
        xtal['2rh1']=value
        file='./prod2/new-matlab-3p0g-%s-aucs-95ci.dat' % type
        value=numpy.loadtxt(file, usecols=(0,))
        xtal['3p0g']=value
        
        best_aucs=[]
        best_paths=[]
        tally=[]
        ohandle=open('./prod2/all-%s-%s-state-aucs.dat' % (sys, type), 'w')
        numpaths=len(glob.glob('./prod2/pathways/sort-%s-path*-states.txt' % (sys)))
        oldmax=0
        for entry in reversed(sorted(best.iteritems(), key=operator.itemgetter(1))):
            for path in numpy.arange(0, numpaths):
                states=numpy.loadtxt('./prod2/pathways/sort-%s-path%s-states.txt' % (sys, path), dtype=int)
                if sys=='car' and state=='53' or state=='63':
                    ohandle.write('%s\t%s\t%s\n' % (state, path, 0))
                else:
                    if int(entry[0]) in states:
                        if entry[0] not in tally:
                            ohandle.write('%s\t%s\t%s\n' % (entry[0], path, entry[1]))
                            tally.append(entry[0])
                        else:
                            pass
                    else:
                        pass
    
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys)
 
