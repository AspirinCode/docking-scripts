#!/bin/python
import numpy
import glob
import optparse
import operator
import random
import os


types=['agonist', 'antagonist']
systems=['2rh1', '3p0g']
for sys in systems:
    for type in types:
        file='./xtal-docking/%s/ranked-%s-%sligand.out' % (type, sys, type)
        ligands=numpy.loadtxt(file, usecols=(1,))
        print "start %s ligands" % (len(ligands))
        ligandnames=numpy.loadtxt(file, usecols=(0,), dtype=str)
        stereo=dict()
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
        file='./xtal-docking/%s/ranked-%s-xtalligand%s-decoys.out' % (type, sys, type)
        decoys=numpy.loadtxt(file, usecols=(1,))
        ofile=open('unsorted.tmp', 'w')
        for x in ligands:
            ofile.write('1\t%s\n' % x)
        for x in decoys:
            ofile.write('0\t%s\n' % x)
        ofile.close()
        os.system('sort -rnk2 unsorted.tmp > ./xtal-docking/%s/new-fmt-%s-%ss.out' % (type, sys, type))
        print "final length %s" % (len(decoys)+len(ligands)) 

