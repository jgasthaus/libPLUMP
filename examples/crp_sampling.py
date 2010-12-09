"""Some utility functions for testing the dynamic programming sampler for 
sampling from CRP_ct."""

import libplump as lp
from collections import defaultdict

def approx_crp_ct_probs_fb(d,c,t,N=100000):
    l = []
    for i in xrange(N):
        l.append(lp.sample_crp_ct_fb(d,c,t))
    d = defaultdict(int)
    for x in l: d[x] += 1
    for x in d: d[x] = d[x]/float(N)
    return d


def approx_crp_ct_probs_bf(d,c,t,N=100000):
    l = []
    for i in xrange(N):
        l.append(lp.sample_crp_ct_bf(d,c,t))
    d = defaultdict(int)
    for x in l: d[x] += 1
    for x in d: d[x] = d[x]/float(N)
    return d

