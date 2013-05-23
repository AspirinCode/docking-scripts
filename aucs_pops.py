import numpy
import pylab
import optparse

def main(sys, type):
    colors=['black', 'grey', 'red', 'orange', 'yellow', 'cyan', 'blue',
            'violet', 'magenta', 'green' ]
    colors=[colors]*4
    pops=numpy.loadtxt('./Populations.dat')
    file='%s_paths.txt' % sys
    fhandle=open(file)
    fluxes=numpy.loadtxt('%s_fluxes.txt' % sys)
    allstates=[]
    pylab.figure()
    for (n, path) in enumerate(fhandle.readlines()):
        states=[int(x) for x in path.split()]
        aucs=[]
        select_states=[]
        for state in states:
            if state not in allstates:
                auc=numpy.loadtxt('%s/new-matlab-bi-%s-%s-aucs-95ci.dat' % (type, state, type), usecols=(0,))
                aucs.append(auc)
                allstates.append(state)
                select_states.append(state)
        if len(select_states)!=0:
            frames=numpy.array(select_states)
            statepops=pops[frames]
            pylab.scatter([-0.6*numpy.log(i) for i in statepops], aucs, s=50, label='path%s' % n, c=colors[n])
            #pylab.scatter(statepops, aucs, c=colors[n])
            pylab.hold(True)
    if type=='agonist':
        xtal='3p0g'
        name='active xtal'
    if type=='antagonist':
        xtal='2rh1'
        name='inactive xtal'
    ref=numpy.loadtxt('../xtal-docking/%s/auc-new-fmt-%s-%s.out' % (type, xtal, type), usecols=(0,))
    pylab.scatter(0.5, ref, marker='x', s=100,  label=name)
    pylab.legend()
    pylab.xlabel('Free Energy')
    pylab.ylim(0.6,1.0)
    pylab.xlim(0,8)
    pylab.ylabel('State AUC')
    pylab.savefig('%s_%s_state_auc_pop.png' % (sys, type))
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

