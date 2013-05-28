from msmbuilder import Trajectory, Project, Serializer
import optparse
import numpy
import os
from numpy import linalg


def main(sys, k, lag, type):
    dir='/Users/morganlawrenz/Desktop/GPCRs/wontkill/new-results/%s/structural-data/' % sys
    types=['path', 'new-path']
    h36=numpy.loadtxt('%s/%s%s.h36.dat' % (dir, type, i))
    npxxy=numpy.loadtxt('%s/inactive-npxxy-%s%s.dat' % (dir, type, i))
    conn=numpy.loadtxt('%s/inactive-conn-%s%s.dat' % (dir, type, i))
    bulge=numpy.loadtxt('%s/active-h5buldge-%s%s.dat' % (dir, type, i))
    dir='%s/Model_l%s' % (dir, lag)
    map=numpy.loadtxt('%s/Mapping.dat' % dir)
    #map=numpy.loadtxt('/home/shukla/GPCR_Mar27/b2ar3p0g2rh1_car/Datak3000/Model_l2/Mapping.dat')
    testactive=numpy.loadtxt('%s/TPT-strict/active_states.txt' % dir)
    testinactive=numpy.loadtxt('%s/TPT-strict/inactive_states.txt' % dir)
    if type=='leu':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0)&(leu75_pro323<8.0))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(leu75_pro323>9.6))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ac_npxxy>1.0)&(ac_conn>1.0)&(leu75_pro323>8.0))[0]
    if type=='asn':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0)&(asn322_ser120<6.0))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(asn322_ser120>8.5))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ac_npxxy>1.0)&(ac_conn>1.0)&(asn322_ser120>8.5))[0]
    if type=='tyr':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0)&(tyr219_tyr326<5.9))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(tyr219_tyr326>11.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ac_npxxy>1.0)&(ac_conn>1.0)&(tyr219_tyr326>11.0))[0]
    if type=='ile':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0)&(ile121_met215>6.0)&(ile121_met215<10.0))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ile121_met215<3.9))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ac_npxxy>1.0)&(ac_conn>1.0)&(ile121_met215<5.0))[0]
    if type=='phe':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0)&(phe208_phe282<5.0))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(phe208_phe282>7.1))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ac_npxxy>1.0)&(ac_conn>1.0)&(phe208_phe282>7.1))[0]
    if type=='arg':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0)&(arg131_glu268>16.4))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0)&(arg131_glu268<7.5))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(ac_npxxy>1.0)&(in_conn<1.0)&(ac_conn>1.0)&(arg131_glu268<7.5))[0]
    if type=='strict-new2':
        frames=numpy.where((h36>12.0)&(h5buldge<1.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(h5buldge>1.0)&(ac_npxxy>1.0)&(in_npxxy<0.8)&(in_conn<1.0)&(ac_conn>1.0))[0]
        #frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        #frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(ac_npxxy>1.0)&(in_conn<1.0)&(ac_conn>1.0))[0]
    elif type=='strict':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(ac_npxxy>1.0)&(in_conn<1.0)&(ac_conn>1.0))[0]
    elif type=='h5':
        frames=numpy.where((h36>12.0)&(h5buldge<1.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(h5buldge>1.0)&(in_npxxy<0.8)&(ac_npxxy>1.0)&(in_conn<1.0)&(ac_conn>1.0))[0]
    elif type=='medium':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(ac_conn<1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0))[0]
        #frames=numpy.where((h36>12.0)&(ac_npxxy<1.0)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<1.0)&(ac_npxxy>1.0)&(in_conn<1.0)&(ac_conn>1.0))[0]
    elif type=='connector-first':
        #frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        frames=numpy.where((h36>12.0)&(in_npxxy>0.8)&(in_conn>1.0))[0]
        inactiveframes=numpy.where((h36>9.0)&(in_npxxy>0.8)&(in_conn<1.0))[0]
    elif type=='connector-second':
        frames=numpy.where((h36>12.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
        inactiveframes=numpy.where((h36>12.0)&(in_npxxy<0.8)&(in_conn<1.0))[0]
    elif type=='npxxy-first':
        frames=numpy.where((h36>9.0)&(in_npxxy>1.0)&(in_conn<1.0))[0]
        inactiveframes=numpy.where((h36>9.0)&(in_npxxy<0.8)&(in_conn<1.0))[0]
    elif type=='npxxy-second':
        frames=numpy.where((h36>9.0)&(in_npxxy>0.8)&(in_conn<1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy>0.8)&(in_conn<1.0))[0]
    elif type=='connector-h36':
        frames=numpy.where((h36>9.0)&(in_npxxy<0.8)&(in_conn>1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0))[0]
    elif type=='npxxy':
        frames=numpy.where((h36>12.0)&(in_npxxy>1.0))[0]
        inactiveframes=numpy.where((h36>12.0)&(in_npxxy<0.8))[0]
    elif type=='h36':
        frames=numpy.where((h36>12.0)&(in_npxxy<0.8)&(in_conn<1.0))[0]
        inactiveframes=numpy.where((h36<9.0)&(in_npxxy<0.8)&(in_conn<1.0))[0]
    #frames=numpy.where((h36>14.0)&(ac_npxxy<0.8)&(in_npxxy>0.8)&(ac_conn<1.0)&(in_conn>1.0))[0]
    #inactiveframes=numpy.where((h36<9.0)&(in_npxxy<1.0)&(ac_npxxy>1.0)&(in_conn<1.0)&(ac_conn>1.0))[0]
    if len(frames)==0:
        import pdb
        pdb.set_trace()
    if len(inactiveframes)==0:
        import pdb
        pdb.set_trace()
    if not os.path.exists('%s/TPT-%s/' % (dir, type)):
        os.mkdir('%s/TPT-%s/' % (dir, type))
    ohandle1=open('%s/TPT-%s/active_states.txt' % (dir, type), 'w')
    ohandle2=open('%s/TPT-%s/gen_active_states.txt' % (dir, type), 'w')
    ihandle1=open('%s/TPT-%s/inactive_states.txt' % (dir, type), 'w')
    ihandle2=open('%s/TPT-%s/gen_inactive_states.txt' % (dir, type), 'w')
    count=0
    for x in frames:
        if map[x] !=-1:
            if map[x] in testactive:
                count+=1
                ohandle1.write('%s\n' % int(map[x]))
                ohandle2.write('%s\n' % int(x))
    print "%s %s active states match with %s in strict" % (type, len(inactiveframes), count)
    count=0
    for x in inactiveframes:
        if map[x] !=-1:
            if map[x] in testinactive:
                count+=1
                ihandle1.write('%s\n' % int(map[x]))
                ihandle2.write('%s\n' % int(x))
    print "%s %s inactive states match with %s in strict" % (type, len(inactiveframes), count)


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-l', '--lag', dest='lag',
                      help='lagtime')
    parser.add_option('-s', '--system', dest='sys',
                      help='system name')
    parser.add_option('-k', '--kclusters', dest='k',
                      help='number of clusters')
    parser.add_option('-t', '--type', dest='type',
                      help='TPT type')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys,k=options.k, lag=options.lag, type=options.type)

