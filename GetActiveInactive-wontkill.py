import optparse
import pylab
import numpy
import os
from numpy import linalg

def op_path(op, name, op_scores, score):
    test=sorted(op_scores[name].keys())
    for (n, x) in enumerate(op):
        for (m, i) in enumerate(test):
            if x >= test[m]:
                if (m+1) >= len(test):
                    score[n]+=op_scores[name][i]
                elif x < test[m+1]:
                    score[n]+=op_scores[name][i]
                    break
    return score

def get_scores():
    max_score=0
    op_scores=dict()
    op_scores['h36']=dict()
    op_scores['conn']=dict()
    op_scores['npxxy']=dict()
    op_scores['bulge']=dict()
    vals=numpy.arange(8, 13,1)
    scores=range(0,len(vals))
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['h36'][v]=s
    vals=numpy.arange(0.5, 3.0, 0.5)
    scores=range(0,len(vals))[::-1]
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['conn'][v]=s
    vals=numpy.arange(0.5, 2.0, 0.5)
    scores=range(0,len(vals))
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['npxxy'][v]=s
    vals=numpy.arange(0.5, 2.5, 0.5)
    scores=range(0,len(vals))[::-1]
    max_score+=max(scores)
    for (v,s) in zip(vals, scores):
        op_scores['bulge'][v]=s
    return max_score, op_scores

def main(sys, aucname):
    dir='/Users/morganlawrenz/Desktop/GPCRs/wontkill/new-results/%s/structural-data/' % sys
    types=['path', 'new-path']
    type='path'
    states=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(0,), dtype=int)
    aucs=numpy.loadtxt('./%s/%s_%s_auc_ci.txt' % (sys, sys, aucname), usecols=(1,))
    paths=open('%s/%s_paths.txt' % (sys, sys))
    map=numpy.loadtxt('%s/Mapping.dat' % dir)
    pop=numpy.loadtxt('%s/Populations.dat' % dir)
    ops=dict()
    all_aucs=dict()
    max_score, op_scores=get_scores()
    print "max score is ", max_score
    all_active=[]
    all_inactive=[]
    for (i, path) in enumerate(paths.readlines()):
        path=numpy.array(path.split())
        score=numpy.zeros(len(path))
        h36=numpy.loadtxt('%s/%s%s.h36.dat' % (dir, type, i))
        bulge=numpy.loadtxt('%s/active-h5buldge-%s%s.dat' % (dir, type, i))
        in_npxxy=numpy.loadtxt('%s/inactive-npxxy-%s%s.dat' % (dir, type, i))
        ac_npxxy=numpy.loadtxt('%s/active-npxxy-%s%s.dat' % (dir, type, i))
        in_conn=numpy.loadtxt('%s/inactive-conn-%s%s.dat' % (dir, type, i))
        ac_conn=numpy.loadtxt('%s/active-conn-%s%s.dat' % (dir, type, i))
        score=op_path(h36, 'h36', op_scores, score)
        score=op_path(ac_conn, 'conn', op_scores, score)
        score=op_path(in_npxxy, 'npxxy', op_scores, score)
        score=op_path(bulge, 'bulge', op_scores, score)
        ops[i]=score
        all_aucs[i]=[]
        for state in path:
            location=numpy.where(states==int(state))[0]
            all_aucs[i].append(aucs[location][0])
    histo=dict()
    for path in sorted(ops.keys()):
        for (n, score) in enumerate(ops[path]):
            if score not in histo.keys():
                histo[score]=[]
                histo[score].append(all_aucs[path][n])
            else:
                histo[score].append(all_aucs[path][n])
    pylab.figure()
    for score in sorted(histo.keys())[::-1]:
        pylab.plot([score]*len(histo[score]), histo[score])
    pylab.show()


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system name')
    parser.add_option('-n', '--aucname', dest='aucname',
                      help='auc type name')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, aucname=options.aucname)


