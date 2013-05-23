#!/bin/python
import numpy
import glob
import optparse
import operator
import random
import os
import operator


def main(sys):
    states=numpy.loadtxt('%s_states.txt' % sys)
    types=['antagonist']
    for type in types:
        n=1
        for state in states:
            if not os.path.exists('%s/%s-state%s/' % (type, sys, n)):
                os.mkdir('%s/%s-state%s/' % (type, sys, n))
            print "on %s for %s" % (state, type)
            file='%s/ranked-%s-aln-state%s.cent-ligand%s.out'% (type, sys,
                    int(state), type) 
            if os.path.exists(file):
                os.system('cp %s %s/%s-state%s/' % (file, type, sys, n))
            file=open('%s/%s-state%s/ranked-%s-aln-state%s.cent-ligand%s.out' %
                    (type, sys, n, sys, int(state), type))
            stereo=dict()
            ligands=[]
            ligandnames=[]
            for line in file.readlines():
                if '::::' not in line:
                    if 'file' not in line:
                        ligands.append(float(line.split()[1]))
                        ligandnames.append(line.split()[0].split('/')[2].split(':')[1])
            ligands=numpy.array(ligands)
            print "start %s ligands" % (len(ligands))
            for (m,x) in enumerate(ligandnames):
                name=x.split('stereo-')[1].split('-')[0]
                if name in stereo.keys():
                    if stereo[name]<float(ligands[m]):
                        stereo[name]=float(ligands[m])
                    else:
                        pass
                else:
                    stereo[name]=float(ligands[m])
            print "final %s ligands" % (len(stereo.keys()))
            file='%s/ranked-%s-aln-state%s.cent-ligand%s-decoys.out'% (type, sys, int(state), type) 
            if os.path.exists(file):
                os.system('cp %s %s/%s-state%s/' % (file, type, sys, n))
            file=open('%s/%s-state%s/ranked-%s-aln-state%s.cent-ligand%s-decoys.out' % (type, sys, n, sys, int(state), type))
            decoys=[]
            for line in file.readlines():
                if '::::' not in line:
                    if '.out' not in line:
                        if 'could not' not in line:
                            decoys.append(float(line.split()[1]))
            decoys=numpy.array(decoys)
            ofile=open('unsorted.tmp', 'w')
            sorted_stereo=reversed(sorted(stereo.iteritems(),
                key=operator.itemgetter(1)))
            for x in sorted_stereo:
                ofile.write('1\t%s\n' % x[1])
            for x in decoys:
                ofile.write('0\t%s\n' % x)
            ofile.close()
            os.system('sort -rnk2 unsorted.tmp > %s/%s-state%s/new-fmt-%s.out' %
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
