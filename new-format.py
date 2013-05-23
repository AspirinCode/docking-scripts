import numpy
import glob
import optparse
import operator
import random
import os
import operator


def main(sys):
    states=numpy.loadtxt('%s_states.txt' % sys)
    types=['agonist' , 'antagonist']
    for type in types:
        n=1
        for state in states:
            all=dict()
            if not os.path.exists('%s/%s-state%s/' % (type, sys, n)):
                os.mkdir('%s/%s-state%s/' % (type, sys, n))
            print "on %s for %s" % (state, type)
            file='%s/ranked-%s-aln-state%s.cent-ligand%s.out'% (type, sys,
                    int(state), type) 
            if os.path.exists(file):
                os.system('cp %s %s/%s-state%s/' % (file, type, sys, n))
            else:
                print "no file"
                import pdb
                pdb.set_trace()
            file=open('%s/%s-state%s/ranked-%s-aln-state%s.cent-ligand%s.out' %
                    (type, sys, n, sys, int(state), type))
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
                all[map[x]]=stereo[x]
            numpy.savetxt('all_ligands.txt', stereo.keys(), fmt='%s')
            file='%s/ranked-%s-aln-state%s.cent-ligand%s-decoys.out'% (type, sys, int(state), type) 
            if os.path.exists(file):
                os.system('cp %s %s/%s-state%s/' % (file, type, sys, n))
            else:
                print "no file"
                import pdb
                pdb.set_trace()
            file=open('%s/%s-state%s/ranked-%s-aln-state%s.cent-ligand%s-decoys.out' % (type, sys, n, sys, int(state), type))
            decoys=dict()
            for line in file.readlines():
                if len(line.split())>1:
                    name=line.split()[0].split('prot-')[1].split('_000')[0]
                    decoys[name]=float(line.split()[1])
                    all[name]=float(line.split()[1])
            ofile=open('unsorted.tmp', 'w')
            sorted_stereo=reversed(sorted(stereo.iteritems(),
                key=operator.itemgetter(1)))
            sorted_decoys=reversed(sorted(decoys.iteritems(),
                key=operator.itemgetter(1)))
            for x in sorted_stereo:
                ofile.write('1\t%s\n' % x[1])
            for x in sorted_decoys:
                ofile.write('0\t%s\n' % x[1])
            ofile.close()
            sorted_all=reversed(sorted(all.iteritems(),
                key=operator.itemgetter(1)))
            num=len(all.keys())*0.2
            selection=[]
            count=0
            for x in sorted_all:
                if count < num:
                    selection.append(x[0])
                    count+=1
                else:
                    break
            numpy.savetxt('top_2percent.txt', selection, fmt='%s')
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
