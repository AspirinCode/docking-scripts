import numpy
import glob
import optparse
import operator
import random
import os
import operator


def main(sys):
    if sys=='apo':
        target=103
    if sys=='bi':
        target=21
    if sys=='car':
        target=44
    states=range(1,target)
    if not os.path.exists('eval-types/'):
        os.mkdir('eval-types/')
    types=['agonist' , 'antagonist']
    n=1
    for state in states:
        all=dict()
        for type in types:
            all[type]=dict()
            if not os.path.exists('eval-types/%s-state%s/' % ( sys, n)):
                os.mkdir('eval-types/%s-state%s/' % (type, sys, n))
            print "on %s for %s" % (state, type)
            file='%s/ranked-%s-aln-state%s.rand-allligand%s.out'% ( type, sys,
                    int(state), type) 
            if os.path.exists(file):
                os.system('cp %s eval-types/%s-state%s/' % (file,  sys, n))
            else:
                print "no file"
                import pdb
                pdb.set_trace()
            file=open('eval-types/%s-state%s/ranked-%s-aln-state%s.rand-allligand%s.out' %
                    ( sys, n, sys, int(state), type))
            stereo=dict()
            map=dict()
            ligands=[]
            ligandnames=[]
            for line in file.readlines():
                if '::::' not in line:
                    if 'file' not in line:
                        ligands.append(float(line.split()[1]))
                        ligandnames.append(line.split()[0])
            ligands=numpy.array(ligands)
            print "start %s ligands" % (len(ligands))
            for (m,x) in enumerate(ligandnames):
                name=x.split('stereo-')[1].split('-')[0]
                if name in stereo.keys():
                    if stereo[name]<float(ligands[m]):
                        stereo[name]=float(ligands[m])
                        map[name]=x.split('_000')[0].split('stereo-')[1]
                    else:
                        pass
                else:
                    stereo[name]=float(ligands[m])
                    map[name]=x.split('_000')[0].split('stereo-')[1]
            print "final %s ligands" % (len(stereo.keys()))
            for x in stereo.keys():
                all[type][map[x]]=stereo[x]
        ofile=open('unsorted.tmp', 'w')
        sorted_agonist=reversed(sorted(all['agonist'].iteritems(),
            key=operator.itemgetter(1)))
        sorted_antagonist=reversed(sorted(all['antagonist'].iteritems(),
            key=operator.itemgetter(1)))
        for x in sorted_agonist:
            ofile.write('1\t%s\n' % x[1])
        for x in sorted_antagonist:
            ofile.write('0\t%s\n' % x[1])
        ofile.close()
        os.system('sort -rnk2 unsorted.tmp > eval-types/%s-state%s/new-fmt-rand-eval-types.out' %
            (sys, n))
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
