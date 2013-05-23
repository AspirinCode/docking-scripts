import numpy
import pylab
import optparse

def main(sys, type):
    colors=['#0000ff','#00ccff','#00ffcc','#ff0000','#ff00cc','#ff66ff','#ff6600','#ff9900','#6600cc','#330099','#00ff00']
    pops=numpy.loadtxt('./Populations.dat')
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    allaucs=[]
    relativef=[]
    allstates=[]
    pylab.figure()
    colorn=0
    for (n, path) in enumerate(fhandle.readlines()):
        states=[int(x) for x in path.split()]
        aucs=[]
        select_states=[]
        flux=fluxes[n]/fluxes[0]
        for state in states:
            if state not in allstates:
                auc=numpy.loadtxt('%s/new-matlab-%s-%s-%s-aucs-95ci.dat' %
                        (type, sys, state, type), usecols=(0,))
                allaucs.append(auc)
                allstates.append(state)
                relativef.append(flux)
    frames=numpy.array(allstates)
    statepops=pops[frames]
    pylab.scatter([-0.6*numpy.log(i) for i in statepops], allaucs,
                        marker='o', s=50, c=relativef,
                        cmap=pylab.get_cmap("hot_r"))
    pylab.colorbar()
    xtal='2rh1'
    name='inactive xtal'
    ref=numpy.loadtxt('../xtal-docking/%s/auc-new-fmt-%s-%s.out' % (type, xtal, type), usecols=(0,))
    refstates=numpy.loadtxt('inactive_xtal_msm_state.txt', dtype=int)
    refpops=pops[refstates]
    pylab.scatter([-0.6*numpy.log(i) for i in refpops], [ref]*len(refstates),
            marker='x', s=100,  label=name)
    xtal='3p0g'
    name='active xtal'
    ref=numpy.loadtxt('../xtal-docking/%s/auc-new-fmt-%s-%s.out' % (type, xtal, type), usecols=(0,))
    refstates=numpy.loadtxt('xtal_msm_state.txt', dtype=int)
    refpops=pops[refstates]
    pylab.scatter([-0.6*numpy.log(i) for i in refpops], [ref]*len(refstates),
            marker='o', s=100,  label=name)
    pylab.legend(loc=2)#bbox_to_anchor=(0.2, 1.1))
    pylab.xlabel('Free Energy')
    pylab.ylim(0.6,1.0)
    pylab.xlim(0,8)
    pylab.ylabel('State AUC')
    pylab.savefig('%s_%s_state_flux_auc_pop.png' % (sys, type))
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    parser.add_option('-t', '--type', dest='type',
            help='type: agonist, antagonist')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, type=options.type)

