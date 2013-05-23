#!/bin/python
import numpy
import glob
import optparse
import operator
import random
import os
import operator


def main(sys):
    #states=numpy.loadtxt('sample_rando_states.txt', dtype=int)
    states=range(1, 21)
    types=['agonist', ] # 'antagonist']
    for type in types:
        n=1
        for state in states:
            if not os.path.exists('%s/%s-state%s/' % (type, sys, n)):
                os.mkdir('%s/%s-state%s/' % (type, sys, n))
            print "on %s for %s" % (state, type)
            file='%s/ranked-%s-aln-state%s.rand-allligand%s.out'% (type, sys,
                    int(state), type) 
            if os.path.exists(file):
                os.system('cp %s %s/%s-state%s/' % (file, type, sys, n))
            else:
                print "no file"
                import pdb
                pdb.set_trace()
            sorted_stereo=numpy.loadtxt(file, usecols=(1,))
            print "final %s ligands" % len(sorted_stereo)
            file='%s/ranked-%s-aln-state%s.rand-ligand%s-decoys.out'% (type, sys, int(state), type) 
            if os.path.exists(file):
                os.system('cp %s %s/%s-state%s/' % (file, type, sys, n))
            else:
                print "no file"
                import pdb
                pdb.set_trace()
            sorted_decoys=numpy.loadtxt(file, usecols=(1,))
            ofile=open('unsorted.tmp', 'w')
            for x in sorted_stereo:
                ofile.write('1\t%s\n' % x)
            for x in sorted_decoys:
                ofile.write('0\t%s\n' % x)
            ofile.close()
            os.system('sort -rnk2 unsorted.tmp > %s/%s-state%s/new-fmt-rand-%s.out' %
                    (type, sys, n, type))
            n+=1

        
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys)
