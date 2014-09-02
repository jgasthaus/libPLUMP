import libplump
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import time
from collections import defaultdict
from antaresia.nputils import * 
from antaresia.plot import * 



log_stirling_cache = {}
def log_stirling_cached(d,n,m):
  d,n,m = (float(d), int(n), int(m))
  if (d,n,m) in log_stirling_cache:
    return log_stirling_cache[(d,n,m)]
  else:
    r = libplump.log_gen_stirling_direct(d,n,m)
    log_stirling_cache[(d,n,m)] = r
    return r

# stratified resampling of Carpenter et al.
def strat_res(w,N):
    K = np.sum(w/N)
    U = np.random.rand()*K
    out = []
    for i in range(len(w)):
      U = U - w[i]
      if U < 0:
          out.append(i)
          U = U + K
    return out


# comp optimal threshold for resampling as in Fearnhead
def opt_thresh(w, N):
    w_sort = np.sort(w)
    o = np.ones((2,len(w)))
    j = -1
    for i in range(len(w)):
      o[1,:] = w_sort/w_sort[i]
      m = np.sum(np.min(o,0))
      if m <= N:
        j = i
        break
    if j == -1:
        print "warning, no such j found!"
    kappa = w_sort[j]
    print kappa
    Ak = np.sum(w_sort >= kappa)
    print Ak
    Bk = np.sum(w_sort[w_sort < kappa])
    return Bk/(N - Ak)

# perform Fearnhead's threshold resampling
def threshold_resampling(w, N):
    outweights = np.zeros(w.shape)
    c = opt_thresh(w,N)
    keep = w > c
    M = np.sum(keep)
    print M
    outweights[keep] = w[keep]
    otherweights = w[w <= c]
    otherweights /= np.sum(otherweights)
    others = strat_res(otherweights,(N-M))
    idx = np.arange(len(w))[w <= c][others]
    outweights[idx] = c
    return outweights



# compute log posterior probability of tw tables up to 
# and additive constant
def log_post_t(a,d,cw,tw,T,p0):
  # see eq. 10 in "Improvements ..." paper
 return (libplump.logKramp(a + d, d, float(T -1)) + 
         log_stirling_cached(d, cw, tw) + 
         tw*math.log(p0))

def log_post_ts(a,d,cws,tws,p0s):
  # see eq. 10 in "Improvements ..." paper
 return (libplump.logKramp(a + d, d, float(np.sum(tws)) - 1) + 
         np.sum([log_stirling_cached(d,cws[i],tws[i]) for i in range(len(cws))]) + 
         np.sum([tws[i]*math.log(p0s[i]) for i in range(len(cws))])
         )

def logcrp(a,d,c,t):
 return (libplump.logKramp(a + d, d, float(t -1)) -
         libplump.logKramp(a + 1, 1, float(c-1)) +
         log_stirling_cached(d, c, t))

def prior_et(a,d,c):
   return (np.exp( libplump.logKramp(a + d, 1, c)
                 - np.log(d) 
                 - libplump.logKramp(a + 1, 1, c - 1))
         - a/d)

def log_post_hierarchical_ts(a,d,cw1,tw1,tw0,p0):
  """Posterior distribution over the number of tables in a two-level hierarchy, 
  where cw1/tw1 are the number of customer at the bottom level, tw0 
  are the tables at the top level and p0 is the base distribution."""
  return (  logcrp(a,d,cw1,tw1) 
          + logcrp(a,d,tw1,tw0)
          + tw0*math.log(p0))

def post_hierarchical(a,d,N,p0):
    res = np.ones((N,N))*-np.inf
    for tw1 in range(N):
        for tw0 in range(tw1+1):
          res[tw1,tw0] = log_post_hierarchical_ts(a,d,N,tw1+1,tw0+1,p0)
    res -= np.max(res)
    res = np.exp(res)
    res /= np.sum(res)
    return res


def post_hierarchical_three(a,d,N,p0):
    res = np.ones((N,N,N))*-np.inf
    for tw2 in range(N):
        for tw1 in range(tw2+1):
            for tw0 in range(tw1+1):
                res[tw2,tw1,tw0] = (  logcrp(a,d,N,tw2 + 1)
                                    + logcrp(a,d,tw2 + 1,tw1 + 1)
                                    + logcrp(a,d,tw1 + 1, tw0 + 1)
                                    + (tw0 + 1)*math.log(p0))
    res -= np.max(res)
    res = np.exp(res)
    res /= np.sum(res)
    return res

def expect_hierarchical_three(a, d, N, p0):
   etw0 = np.sum(np.arange(1,N+1)*np.sum(np.sum(post_hierarchical_three(a, d, N, p0), 0), 0))
   etw1 = np.sum(np.arange(1,N+1)*np.sum(np.sum(post_hierarchical_three(a, d, N, p0), 2), 0))
   etw2 = np.sum(np.arange(1,N+1)*np.sum(np.sum(post_hierarchical_three(a, d, N, p0), 2), 1))
   return etw2, etw1, etw0

def max_hierarchical_three(a, d, N, p0):
   post = post_hierarchical_three(a, d, N, p0)
   global_idx = np.unravel_index(np.argmax(post), post.shape)
   m0 = np.argmax(np.sum(np.sum(post, 0), 0)) + 1
   m1 = np.argmax(np.sum(np.sum(post, 2), 0)) + 1
   m2 = np.argmax(np.sum(np.sum(post, 2), 1)) + 1
   
   c2 = np.argmax(post[:, global_idx[1], global_idx[2]])
   c1 = np.argmax(post[c2, :, global_idx[2]])
   c0 = np.argmax(post[c2, c1, :])


   return (m2, m1, m0), np.array(global_idx) + 1, (c2 + 1, c1 + 1, c0 + 1)

  
def joint_hierarchical(a,d,N,p0):
    res = np.ones((N,N))*-np.inf
    for tw1 in range(N):
        for tw0 in range(tw1+1):
          res[tw1,tw0] = (logcrp(a,d,N,tw1+1) + logcrp(a,d,tw1+1,tw0+1))
    res = np.exp(res)
    return res

def hierarchical_pf(a,d,N,p0,num_particles=1000, resample=False):
  X = np.zeros((2,num_particles, N))
  W = np.zeros((num_particles, N))
  X[:,:,0] = 1 # first customer at first table
  W[:,0] = 1./num_particles # first customer at first table
  et = np.zeros(N)
  et[0] = 1
  for i in range(1,N):
    for j in range(num_particles):
      cw1 = i
      tw1 = X[0,j,i-1]
      tw0 = X[1,j,i-1]
      p00 = (tw1 - tw0*d + (a+tw0*d)*p0)/(tw1 + a) # predictive in top rest
      p = (cw1 - tw1*d + (a+tw1*d)*p00)/(cw1 + a) # predictive in bottom rest
      f1 = (a + d * tw1)*p00/((a + d*tw1)*p00 + cw1 - tw1*d) # new table prob bottom
      f0 = (a + d * tw0)*p0/((a + d*tw0)*p0 + tw1 - tw0*d) # new table prob top
      bottom_new_table = int(np.random.rand() < f1)
      top_new_table = 0
      if bottom_new_table:
        top_new_table = int(np.random.rand() < f0)
      X[0,j,i] = tw1 + bottom_new_table
      X[1,j,i] = tw0 + top_new_table
      W[j,i] = W[j,i-1] * p 
    W[:,i] /= np.sum(W[:,i])
    if resample and 1./np.sum(np.square(W[:,i])) < 0.1*num_particles:
      #idx = choice(W[:,i], num_particles)
      idx = stratified_resampling(W[:,i])
      idx = residual_resampling(W[:,i])
      X[:,:,:] = X[:,idx,:]
      W[:,i] = 1./num_particles
    #et[i] = np.mean(np.dot(W[:,i], X[1,:,i]))
  #return np.mean(np.dot(W[:,N-1], X[1,:,N-1])), et
  return X, W

def evaluate_pf_result(X,W,at):
   idx = np.sum(X[:,:,-1] == np.array(at)[:,None],0)==2
   return np.sum(W[idx,-1])



def post_t(a, d, cw, other_T, p0):
  """Compute the posterior probability distribution over the
  number of tables in a two parameter CRP, assuming you have observed
  cw customers of a single type."""
  r = np.array([log_post_t(a,d,cw,t, other_T + t, p0) for t in range(1,cw+1)])
  r -= np.max(r)
  r = np.exp(r)
  r /= np.sum(r)
  return r

def pypPredictive(a,d,cwk,twk,cw,tw,p0):
  return (cwk - twk*d + (a+tw*d)*p0)/(cw + a)

def expected_num_tables_two(a,d,cws,p0s):
  ps = np.zeros((cws[0], cws[1]))
  # compute post. probs for each t
  for i in range(cws[0]):
    for j in range(cws[1]):
      ps[i,j] = log_post_ts(a,d,cws,(i+1, j+1),p0s)
  
  # normalize
  ps = ps - np.max(ps)
  ps = np.exp(ps)/np.sum(np.exp(ps))
  marg_0 = np.sum(ps,1)
  marg_1 = np.sum(ps,0)
  return ps, (np.dot(range(1,cws[0]+1),marg_0), np.dot(range(1,cws[1]+1), marg_1))

def expected_num_tables(a,d,N,p0):
  tables = range(1,N+1)
  # compute post. probs for each t
  ps = [log_post_t(a,d,N,t,t,p0) for t in tables]
  
  # normalize
  ps = np.array(ps)
  ps = ps - np.max(ps)
  ps = np.exp(ps)/np.sum(np.exp(ps))
  
  return (ps, np.dot(ps, tables))

def choice(p, N):
  u = np.random.rand(N)
  c = np.cumsum(p)
  return c.searchsorted(u * c[-1])

def optimal_proposal(a, d, cwk, twk, tw, p0):
  post = post_t(a,d,cwk,tw - twk, p0)
  sample = choice(post,1)[0]
  return (sample + 1, post[sample])

def collapsed_particle_num_tables(a, d, obs, p0s, 
    proposal=lambda a,d,cwk,twk,tw,p0: (1 + np.random.randint(cwk), 1./cwk), # uniform proposal
    num_particles=1000, resample=False):
  X = np.zeros((max(obs)+1,2,num_particles),dtype=np.int)
  W = np.zeros(num_particles)
  #X[obs[0],0,:] = 1 # first customer at first table
  #
  #X[obs[0],1,:] = 1 # first customer at first table
  W[:] = 1./num_particles # first customer at first table
  et = np.zeros(len(obs))
  et[0] = 1
  for i in range(len(obs)):
    for j in range(num_particles):
      k = obs[i]
      cwk = X[k,0,j]
      twk = X[k,1,j]
      tw = np.sum(X[:,1,j])
      (tw_proposal, proposal_prob) = proposal(a,d,cwk + 1, twk, tw, p0s[k])
      assert(tw_proposal > 0 and tw_proposal <= cwk+1)
      #print tw_proposal, proposal_prob
      X[k,0,j] += 1 # seat customer 
      X[k,1,j] = tw_proposal
      joint_new = np.exp(log_post_t(a,d, cwk + 1, tw_proposal, tw - twk + tw_proposal, p0s[k]))
      joint_old = np.exp(log_post_t(a,d, cwk, twk, tw, p0s[k]))
      W[j] *= (joint_new/(proposal_prob*joint_old))
    #if resample and 1./np.sum(np.square(W[:,i])) < 20:
    #  W[:,i] /= np.sum(W[:,N-1])
    #  idx = choice(W[:,i], num_particles)
    #  X[:,:,:] = X[:,idx,:]
    #  W[:,i] = 1./num_particles
    W[:] /= np.sum(W[:])
    et[i] = np.dot(W, X[0,1,:])
  return np.dot(W[:], X[:,1,:].T), et


# collapsed particle filter operating on the space of table counts
# this is effectively equivalent to doing straight importance sampling
# for the target distribution
def collapsed_particle_num_tables_single(a, d, N, p0, 
    proposal=lambda cw,tw: (1 + np.random.randint(cw), 1./cw), # uniform proposal
    num_particles=1000, resample=False):
  X = np.zeros((2,num_particles, N))
  W = np.zeros((num_particles, N))
  X[:,:,0] = 1 # first customer at first table
  W[:,0] = 1./num_particles # first customer at first table
  et = np.zeros(N)
  et[0] = 1
  for i in range(1,N):
    for j in range(num_particles):
      cw = X[0,j,i-1]
      tw = X[1,j,i-1]
      (tw_proposal, proposal_prob) = proposal(cw + 1, tw)
      #print tw_proposal, proposal_prob
      X[0,j,i] = X[0,j,i-1] + 1 
      X[1,j,i] = tw_proposal
      joint_new = np.exp(log_post_t(a,d, cw + 1, tw_proposal, tw_proposal, p0))
      joint_old = np.exp(log_post_t(a,d, cw, tw, tw, p0))
      W[j,i] = W[j,i-1] * (joint_new/(proposal_prob*joint_old))
    W[:,i] /= np.sum(W[:,i])
    if resample and 1./np.sum(np.square(W[:,i])) < 20:
      idx = choice(W[:,i], num_particles)
      X[:,:,:] = X[:,idx,:]
      W[:,i] = 1./num_particles
    et[i] = np.mean(np.dot(W[:,i], X[1,:,i]))
  W[:,N-1] /= np.sum(W[:,N-1])
  return np.mean(np.dot(W[:,N-1], X[1,:,N-1])), et

def crp_particle_num_tables(a, d, N, p0, num_particles=1000, resample=False):
  X = np.zeros((2,num_particles, N))
  W = np.zeros((num_particles, N))
  X[:,:,0] = 1 # first customer at first table
  W[:,0] = 1./num_particles # first customer at first table
  et = np.zeros(N)
  et[0] = 1
  for i in range(1,N):
    for j in range(num_particles):
      cw = X[0,j,i-1]
      tw = X[1,j,i-1]
      f = (a + d * tw)*p0/((a + d*tw)*p0 + cw - tw*d)
      p = (cw - tw*d + (a+tw*d)*p0)/(cw + a)
      X[0,j,i] = X[0,j,i-1] + 1 
      X[1,j,i] = X[1,j,i-1] + int(np.random.rand() < f)
      W[j,i] = W[j,i-1] * p
    W[:,i] /= np.sum(W[:,i])
    if resample and 1./np.sum(np.square(W[:,i])) < 0.1*num_particles:
      #idx = choice(W[:,i], num_particles)
      idx = stratified_resampling(W[:,i])
      idx = residual_resampling(W[:,i])
      X[:,:,:] = X[:,idx,:]
      W[:,i] = 1./num_particles
    et[i] = np.mean(np.dot(W[:,i], X[1,:,i]))
  return np.mean(np.dot(W[:,N-1], X[1,:,N-1])), et

def full_crp_particle_num_tables(a, d, N, p0, num_particles=1000, resample=False):
  X = np.zeros((N,num_particles, N)) # table_id, particles, input pos
  W = np.zeros((num_particles, N))
  X[0,:,0] = 1 # first customer at first table
  W[:,0] = 1./num_particles # first customer at first table
  et = np.zeros(N)
  et[0] = 1
  for i in range(1,N):
    X[:,:,i] = X[:,:,i-1]
    for j in range(num_particles):
      tables = X[:,j,i]
      num_tables = np.sum(tables > 0)
      probs = np.zeros(num_tables + 1)
      probs[0:num_tables] = tables[:num_tables] - d
      probs[-1] = (a + d*num_tables)*p0
      sample = choice(probs, 1)[0]
      p = (np.sum(tables) - num_tables*d + (a+num_tables*d)*p0)/(np.sum(tables) + a)
      X[sample,j,i] += 1
      W[j,i] = W[j,i-1] * p
    W[:,i] /= np.sum(W[:,i])
    if resample and 1./np.sum(np.square(W[:,i])) < 0.1*num_particles:
      #idx = choice(W[:,i], num_particles)
      idx = stratified_resampling(W[:,i])
      idx = residual_resampling(W[:,i])
      X[:,:,:] = X[:,idx,:]
      W[:,i] = 1./num_particles
    et[i] = np.mean(np.dot(W[:,i], np.sum(X[:,:,i]>0,0)))
  return np.mean(np.dot(W[:,N-1], np.sum(X[:,:,N-1]>0,0))), et


def crp_particle_num_tables_uniform_proposal(a, d, N, p0, num_particles=1000, resample=False):
  X = np.zeros((2,num_particles, N))
  W = np.zeros((num_particles, N))
  X[:,:,0] = 1 # first customer at first table
  W[:,0] = 1./num_particles # first customer at first table
  et = np.zeros(N)
  et[0] = 1
  for i in range(1,N):
    for j in range(num_particles):
      cw = X[0,j,i-1]
      tw = X[1,j,i-1]
      new_table = int(np.random.rand() < 0.5)

      X[0,j,i] = X[0,j,i-1] + 1 
      X[1,j,i] = X[1,j,i-1] + new_table
      if new_table:
        inc_weight = p0 * (a + d*tw)
      else:
        inc_weight = cw - tw*d
      W[j,i] = W[j,i-1] * inc_weight/0.5
    W[:,i] /= np.sum(W[:,i])
    if resample: # and 1./np.sum(np.square(W[:,i])) < 1000:
      idx = choice(W[:,i], num_particles)
      X[:,:,:] = X[:,idx,:]
      W[:,i] = 1./num_particles
      print X, W
    et[i] = np.mean(np.dot(W[:,i], X[1,:,i]))
  return np.mean(np.dot(W[:,N-1], X[1,:,N-1])), et

def crp_particle_num_tables_enumerate(a, d, N, p0, num_particles=1000, merge = True, sample=True):
  particles = [1]
  weights = [1.]
  et = [1]
  for i in range(1,N):
    new_particles = []
    new_weights = []
    for j in range(len(particles)):
      cw = i
      tw = particles[j]
      f = (a + d * tw)*p0/((a + d*tw)*p0 + cw - tw*d)
      p = (cw - tw*d + (a+tw*d)*p0)/(cw + a)
      w = weights[j]
      new_particles.extend([tw, tw+1])
      new_weights.extend([w*p*(1-f), w*p*f])
      #new_weights.extend([w, w*p0])
    weights = np.array(new_weights)
    weights /= np.sum(weights)
    if merge:
      particles=np.array(new_particles)
      u = np.unique(particles)
      weights = np.bincount(particles, weights)[u]
      particles=u
    else:
      particles = new_particles
    #print particles, weights
    if len(particles) > num_particles:
      if sample:
        newweights = threshold_resampling(weights, num_particles)
        idx = newweights > 0
        print np.sum(idx)
        particles = particles[idx]
        weights = newweights[idx]
        print particles.shape, weights.shape
      else:
        idx = np.argsort(weights)[-num_particles:]
        weights = weights[idx]
        weights /= np.sum(weights)
        particles = [particles[k] for k in idx]
    et.append(np.dot(weights, particles))
  particles = np.array(particles)
  return np.dot(weights, particles), et


def crp_fractional_num_tables(a, d, N, p0):
  cw = 0
  tw = 0
  for i in range(N):
    f = (a + d * tw)*p0/((a + d*tw)*p0 + cw - tw*d)
    tw += f
    cw += 1
  return tw


def crp_fractional_hierarchical_two(a, d, N, p0):
  cw1 = 0 # bottom
  tw1 = 0
  tw0 = 0 # top
  for i in range(N):
    p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
    f1 = (a + d * tw1)*p1/((a + d*tw1)*p1 + cw1 - tw1*d)
    f0 = (a + d * tw0)*p0/((a + d*tw0)*p0 + tw1 - tw0*d) * f1
    tw1 += f1
    tw0 += f0
    cw1 += 1
  return tw1, tw0


def crp_fractional_hierarchical_three(a, d, N, p0):
  """Fractional approximate inference in a three-level
  hierarchy."""
  cw2 = 0 # bottom
  tw2 = 0
  tw1 = 0 # middle
  tw0 = 0 # top
  for i in range(N):
    p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
    p2 = (tw2 - d*tw1 + (a + d*tw1)*p1)/(tw2 + a)
    f2 = (a + d * tw2)*p2/((a + d*tw2)*p2 + cw2 - tw2*d)
    f1 = (a + d * tw1)*p1/((a + d*tw1)*p1 + tw2 - tw1*d) * f2
    f0 = (a + d * tw0)*p0/((a + d*tw0)*p0 + tw1 - tw0*d) * f1
    tw2 += f2
    tw1 += f1
    tw0 += f0
    cw2 += 1
  return tw2, tw1, tw0
  
def crp_optimize_hierarchical_three(a, d, N, p0):
  """Iterative local optimization"""
  cw2 = 0 # bottom
  tw2 = 0
  tw1 = 0 # middle
  tw0 = 0 # top
  for i in range(N):
    cw2 += 1
    for j in range(2):
      p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
      p2 = (tw2 - d*tw1 + (a + d*tw1)*p1)/(tw2 + a)
      post2 = post_t(a, d, cw2, 0, p2)
      tw2 = np.argmax(post2) + 1
      post1 = post_t(a, d, tw2, 0, p1)
      tw1 = np.argmax(post1) + 1
      post0 = post_t(a,d, tw1, 0, p0)
      tw0 = np.argmax(post0) + 1
  return tw2, tw1, tw0

def crp_optimize_hierarchical_alt_three(a, d, N, p0):
  """Iterative local optimization"""
  cw2 = 0 # bottom
  tw2 = 0
  tw1 = 0 # middle
  tw0 = 0 # top
  for i in range(N):
    cw2 += 1
    for j in range(10):
      p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
      p2 = (tw2 - d*tw1 + (a + d*tw1)*p1)/(tw2 + a)
      post2 = post_t(a, d, cw2, 0, p2)
      tw2 = np.argmax(post2) + 1
      tw1 = tw2
      tw0 = tw2
  return tw2, tw1, tw0

def crp_optimize_hierarchical_pairs_three(a, d, N, p0):
  """Iterative local optimization"""
  cw2 = 0 # bottom
  tw2 = 0
  tw1 = 0 # middle
  tw0 = 0 # top
  for i in range(N):
    cw2 += 1
    for j in range(1):
      p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
      p2 = (tw2 - d*tw1 + (a + d*tw1)*p1)/(tw2 + a)
      post1 = post_hierarchical(a, d, cw2, p1)
      tw2, tw1 = np.array(np.unravel_index(np.argmax(post1), post1.shape)) + 1
      post2 = post_hierarchical(a, d, tw2, p0)
      tw1, tw0 = np.array(np.unravel_index(np.argmax(post2), post2.shape)) + 1

  return tw2, tw1, tw0


def post_direct(a, d, cw1, tw0):
  out = np.zeros(cw1)
  for i in range(cw1):
    out[i] = (  libplump.logKramp(a + d, d, i)
              - libplump.logKramp(a + 1, 1, i)
              + log_stirling_cached(d, cw1, i + 1)
              + log_stirling_cached(d, i + 1, tw0))
  out -= np.max(out)
  out = np.exp(out)
  out /= np.sum(out)
  return out

  

def crp_sample_hierarchical_direct(a, d, N, p0, num_samples, initPf=True):
  """Direct gibbs sampler for two layer hierarchy"""
  cw1 = 0 # bottom
  tw1 = 0 # bottom
  tw0 = 0 # top
  res = []
  if initPf:
    for i in range(N):
      p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
      f1 = (a + d * tw1)*p1/((a + d*tw1)*p1 + cw1 - tw1*d)
      f0 = (a + d * tw0)*p0/((a + d*tw0)*p0 + tw1 - tw0*d)
      if (np.random.rand() < f1):
        tw1 +=1
        if (np.random.rand() < f0):
          tw0 +=1
      cw1 += 1
  else:
    for i in range(N):
      cw1 += 1
      p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
      post = post_direct(a, d, cw1, tw0)
      tw1_new = (choice(post, 1) + 1)[0]
      post = post_t(a, d, tw1_new, 0, p0)
      tw0 = (choice(post, 1) + 1)[0]
      tw1 = tw1_new
  for j in range(num_samples):
    p1 = (tw1 - d*tw0 + (a + d*tw0)*p0)/(tw1 + a)
    p2 = (cw1 - d*tw1 + (a + d*tw1)*p1)/(cw1 + a)
    res.append((p2, tw0, tw1))
    post = post_direct(a, d, cw1, tw0)
    tw1 = (choice(post, 1) + 1)[0]
    post = post_t(a, d, tw1, 0, p0)
    tw0 = (choice(post, 1) + 1)[0]
  return res

def crp_fractional_num_tables_p0s(a, d, customers, p0dict):
  cw = defaultdict(lambda: 0)
  tw = defaultdict(lambda: 0)
  c = 0
  t = 0
  for c in customers:
    f = (a + d * t)*p0dict[c]/((a + d*t)*p0dict[c] + cw[c] - tw[c]*d)
    tw[c] += f
    t += f
    cw[c] += 1
    c += 1
  return tw



def plot_posterior_tables_vary_d(a,ds,N,p0):
  i = 1
  fig = plt.figure()
  for d in ds:
    post, et = expected_num_tables(a,float(d),N,p0)
    ax = fig.add_subplot(len(ds)/3, 3, i)
    ax.bar(range(1,N+1), post, color = "black")
    ax.axis((1,N,0,0.5))
    ax.grid();
    ax.set_title("$d = %.1f$, $E[t_w] = %.2f$" %(d, et))
    i += 1
  fig.suptitle("Posterior distribution of $t_s$; $c_s = %d$, $\\alpha=%.1f$, $H(s)=%.1f$" % (N,a, p0))
  fig.set_size_inches(8,9)
  return fig, ax


def plot_posterior_tables_vary_p(a,d,N,p0s):
  i = 1
  fig = plt.figure()
  for p0 in p0s:
    post, et = expected_num_tables(a,float(d),N,p0)
    ax = fig.add_subplot(len(p0s)/3, 3, i)
    ax.bar(range(1,N+1), post, color = "black")
    ax.axis((1,N,0,0.8))
    ax.grid();
    ax.set_title("$H(s) = %.1f$, $E[t_w] = %.2f$" %(p0, et))
    i += 1
  fig.suptitle("Posterior distribution of $t_s$; $c_s = %d$, $\\alpha=%.1f$, $d=%.1f$" % (N,a, d))
  fig.set_size_inches(8,9)
  return fig, ax


def plot_posterior_tables_ds_ps(a,N,ds=np.arange(0.0,1,0.1),p0s=np.arange(0.2,1.1,0.1)):
  i = 1
  cm = plt.cm.spectral
  for p0 in p0s:
    plt.subplot(len(ds)/3, 3, i)
    plt.gca().set_color_cycle([cm(k) for k in np.linspace(0, 1, 11)])
    for d in ds:
      post, et = expected_num_tables(a,float(d),N,float(p0))
      plt.plot(range(1,N+1), post)
    plt.gca().set_color_cycle([cm(k) for k in np.linspace(0, 1, 11)])
    for d in ds:
      post, et = expected_num_tables(a,float(d),N,float(p0))
      plt.plot([et],[d],'x')
    plt.axis((1,N,0,1))
    plt.grid();
    plt.title("H(s) = %.2f" %(p0,))
    i += 1
  plt.gcf().suptitle("Posterior distribution of $t_s$; $c_s = %d$, $\\alpha=%.1f$" % (N,a))



def plot_customers_vs_expected_tables_discounts(a, p0, max_customers, ds=np.arange(0.1,1,0.2), legend_loc=2, plot_frac=True, plot_particles=0, particle_filter=crp_particle_num_tables_enumerate):
  customers = range(1,max_customers + 1)
  et_exact = [[expected_num_tables(a,d,c,p0)[1] for c in customers] for d in ds]
  et_exact = np.array(et_exact).T
  et_approx = [[crp_fractional_num_tables(a,d,c,p0) for c in customers] for d in ds]
  et_approx = np.array(et_approx).T
  #plt.plot(customers, et_exact, customers, et_approx, "--")
  plt.plot(customers, et_exact)
  if plot_frac:
    plt.gca().set_color_cycle(mpl.rcParams['axes.color_cycle'])
    plt.plot(customers, et_approx, '--')
  if plot_particles > 0:
    plt.gca().set_color_cycle(mpl.rcParams['axes.color_cycle'])
    #plt.plot(customers, np.array([[crp_particle_num_tables_enumerate(a,d,c,p0,plot_particles,False) for c in customers] for d in ds]).T, '-.')
    plt.plot(customers, np.array([[particle_filter(a,d,c,p0,num_particles=plot_particles) for c in customers] for d in ds]).T, '-.')
  plt.legend(["$d=%.1f$"%(d,) for d in ds], loc=legend_loc)
  plt.xlabel("$c_w$")
  plt.ylabel("$E[t_w]$")
  plt.title("Posterior $E[t_w]$, $\\alpha=%.1f$, $p_0=%.2f$" % (a, p0))
  plt.grid()

def plot_tables_posterior_two(a, d, cs, p0s):
  plt.imshow(expected_num_tables_two(a,d,cs,p0s)[0], extent=(1,cs[0],1,cs[1]), 
             interpolation="none", origin="lower")

def plot_tables_posterior_two_grid(a,d,cs):
  ps = np.arange(0.1,1,0.1)
  for i in range(len(ps)):
    plt.subplot(3,3,i+1)
    plot_tables_posterior_two(a,d,cs, [ps[i], 1-ps[i]])
    plt.title("p = %f"%(ps[i],))


def plot_posterior_components(a,N, ds=np.arange(0.1,1,0.2), ps=np.arange(0.1,1,0.2)):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  # plot stirling numbers
  ps = []
  for d in ds:
    p1, = ax.plot(range(1,N+1), [log_stirling_cached(d,N,i) for i in range(1,N+1)], ':')

  # plot Kramp's symbol
  ax.set_color_cycle(mpl.rcParams['axes.color_cycle'])
  for d in ds:
    p2, = ax.plot(range(1,N+1), [libplump.logKramp(a + d, d, i - 1) for i in range(1,N+1)],'--')
  
  # plot sum of the two
  ax.set_color_cycle(mpl.rcParams['axes.color_cycle'])
  for d in ds:
    p3, = ax.plot(range(1,N+1), [libplump.logKramp(a + d, d, i - 1) + log_stirling_cached(d,N,i) for i in range(1,N+1)])
    ps.append(p3)
  
  #plt.gca().set_color_cycle(mpl.rcParams['axes.color_cycle'])
  #for p in ps:
  #  plt.plot(range(1,N+1), [np.log(p)*i for i in range(1,N+1)],'.')
  
  
  #plt.gca().set_color_cycle(mpl.rcParams['axes.color_cycle'])
  # for p in ps:
  #   plt.plot(range(1,N+1), [np.log(p)*i for i in range(1,N+1)],'.')
  ax.set_xlabel("Number of tables")
  ax.set_title("Log-Posterior Parts, $\\alpha = %.1f$"%(a,))
  l1 = ax.legend([p1, p2, p3], 
      [ "$\log  S_d(c, t)$",
        "$\log [\\alpha + d]_d^{t-1}$",
        "$\log  S_d(c, t) + \log  [\\alpha + d]_d^{t-1}$"
      ], loc = 6)
  ax.legend(ps, ["$d = %.1f$"%(d,) for d in ds], loc = 8, 
      ncol=len(ds) , borderaxespad=0.)
  ax.add_artist(l1)
  ax.grid(True)
  return fig, ax

def get_errorbars(fun, K=10):
    runs = []
    for i in range(K):
      runs.append(fun())
    data = np.array(runs)
    m = np.mean(data,0)
    sd = np.std(data,0)
    return m, sd


def plot_crp_fractional(a, d, N, p0):
   plt.plot(range(1,N+1),[expected_num_tables(a, d,i, p0)[1] for i in range(1, N + 1)],'k',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'b--',linewidth=2)
   plt.grid()
   plt.title(r"Posterior Expected Number of Tables, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(a,d,p0))
   plt.xlabel("# customers of type $s$")
   plt.ylabel('$E[t_s]$')
   plt.xlim(1,N)


def plot_crp_fractional_hier_two(a, d, N, p0):
   plt.plot(range(1,N+1),[np.sum(np.arange(1,i+1)*np.sum(post_hierarchical(a, d, i, p0), 0)) for i in range(1, N + 1)],'k',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_hierarchical_two(a,d,c,p0)[1] for c in range(1,N+1)],'k--',linewidth=2)
   plt.plot(range(1,N+1),[np.sum(np.arange(1,i+1)*np.sum(post_hierarchical(a, d, i, p0), 1)) for i in range(1, N + 1)],'r',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_hierarchical_two(a,d,c,p0)[0] for c in range(1,N+1)],'r--',linewidth=2)
   plt.grid()
   plt.title(r"Posterior Expected Number of Tables, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(a,d,p0))
   plt.xlabel("# customers of type $s$")
   plt.ylabel('$E[t_s]$')
   plt.xlim(1,N)


def plot_crp_fractional_hier_three(a, d, N, p0):
   plt.plot(range(1,N+1),[np.sum(np.arange(1,i+1)*np.sum(np.sum(post_hierarchical_three(a, d, i, p0), 0), 0)) for i in range(1, N + 1)],'k',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_hierarchical_three(a,d,c,p0)[2] for c in range(1,N+1)],'k--',linewidth=2, label="_nolegend_")


   plt.plot(range(1,N+1),[np.sum(np.arange(1,i+1)*np.sum(np.sum(post_hierarchical_three(a, d, i, p0), 2), 0)) for i in range(1, N + 1)],'r',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_hierarchical_three(a,d,c,p0)[1] for c in range(1,N+1)],'r--',linewidth=2, label="_nolegend_")
   
   plt.plot(range(1,N+1),[np.sum(np.arange(1,i+1)*np.sum(np.sum(post_hierarchical_three(a, d, i, p0), 2), 1)) for i in range(1, N + 1)],'b',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_hierarchical_three(a,d,c,p0)[0] for c in range(1,N+1)],'b--',linewidth=2, label="_nolegend_")

   plt.grid()
   plt.title(r"Marginal Posterior $\mathrm{E}[t_s]$, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(a,d,p0))
   plt.xlabel("# customers of type $s$")
   plt.ylabel('$E[t_s]$')
   plt.xlim(1,N)
   plt.legend(["top", "middle", "bottom"], "upper left")


def plot_crp_fractional_p0s(a, d, N):
   p0 = 0.1
   plt.plot(range(1,N+1),[expected_num_tables(a, d,i, p0)[1] for i in range(1, N + 1)],'k',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'k--',linewidth=2, label="_nolegend_")
   p0 = 0.4
   plt.plot(range(1,N+1),[expected_num_tables(a, d,i, p0)[1] for i in range(1, N + 1)],'r',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'r--',linewidth=2, label="_nolegend_")
   p0 = 0.7
   plt.plot(range(1,N+1),[expected_num_tables(a, d,i, p0)[1] for i in range(1, N + 1)],'b',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'b--',linewidth=2, label="_nolegend_")
   p0 = 0.95
   plt.plot(range(1,N+1),[expected_num_tables(a, d,i, p0)[1] for i in range(1, N + 1)],'g',linewidth=2)
   plt.plot(range(1,N+1),[crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'g--',linewidth=2, label="_nolegend_")


   plt.grid()
   plt.title(r"Posterior Expected Number of Tables, $\alpha=%.1f$, $d=%.1f$"%(a,d))
   plt.xlabel("# customers of type $s$")
   plt.ylabel('$E[t_s]$')
   plt.xlim(1,N)
   plt.legend(["H(s) = 0.1",
               "H(s) = 0.4",
               "H(s) = 0.7",
               "H(s) = 0.95"], "upper left")


def plot_maximization(a, d, N, p0):
  m2s = []
  m1s = []
  m0s = []
  ml2s = []
  ml1s = []
  ml0s = []
  ma2s = []
  ma1s = []
  ma0s = []
  g2s = []
  g1s = []
  g0s = []
  for i in range(1, N+1):
    ml2, ml1, ml0 = crp_optimize_hierarchical_three(a, d, i, p0)
    m2, m1, m0 = crp_optimize_hierarchical_pairs_three(a, d, i, p0)
    m2s.append(m2)
    m1s.append(m1)
    m0s.append(m0)
    ml2s.append(ml2)
    ml1s.append(ml1)
    ml0s.append(ml0)
    margmax, globmax, _ = max_hierarchical_three(a, d, i, p0)
    ma2s.append(margmax[0])
    ma1s.append(margmax[1])
    ma0s.append(margmax[2])
    g2s.append(globmax[0])
    g1s.append(globmax[1])
    g0s.append(globmax[2])
  print g0s
  plt.plot(range(1, N + 1), m2s, "r", alpha=0.7, linewidth=2)
  plt.plot(range(1, N + 1), m1s, "r--", label="_nolegend_", linewidth=2, alpha=0.7)
  plt.plot(range(1, N + 1), m0s, "r-.", label="_nolegend_", linewidth=2, alpha=0.7)
  plt.plot(range(1, N + 1), ml2s, "g", alpha=0.7, linewidth=2)
  plt.plot(range(1, N + 1), ml1s, "g--", label="_nolegend_", linewidth=2, alpha=0.7)
  plt.plot(range(1, N + 1), ml0s, "g-.", label="_nolegend_", linewidth=2, alpha=0.7)
  # plt.plot(range(1, N + 1), ma2s,"b", linewidth=2)
  # plt.plot(range(1, N + 1), ma1s,"b--", label="_nolegend_", linewidth=2)
  # plt.plot(range(1, N + 1), ma0s,"b-.", label="_nolegend_", linewidth=2)
  plt.plot(range(1, N + 1), g2s, "k", linewidth=2, alpha=0.7)
  plt.plot(range(1, N + 1), g1s, "k--", label="_nolegend_", linewidth=2, alpha=0.7)
  plt.plot(range(1, N + 1), g0s, "k-.", label="_nolegend_", linewidth=2, alpha=0.7)
  plt.grid()
  plt.title(r"Posterior Number of Tables, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(a,d,p0))
  plt.xlabel("# customers of type $s$")
  plt.ylabel('$t_s$')
  plt.xlim(1,N)
  plt.legend(["pairwise", "local", "global"], "upper left")


def plot_crp_particle_filters(a, d, N, p0, K=10):
   #m1, s1 = get_errorbars(lambda:  full_crp_particle_num_tables(a,d,N,p0,100)[1],K)
   #plt.errorbar(np.arange(1,N+1)-0.3, m1, s1,fmt='b',alpha=0.5)
   m1, s1 = get_errorbars(lambda: crp_particle_num_tables(a,d,N,p0,1)[1], K)
   plt.errorbar(np.arange(1,N+1), m1, s1,fmt='b',alpha=0.5)
   m2, s2 = get_errorbars(lambda: crp_particle_num_tables(a,d,N,p0,100)[1], K)
   plt.errorbar(np.arange(1,N+1), m2, s2,fmt='r',alpha=0.5)
   m3, s3 = get_errorbars(lambda: crp_particle_num_tables_enumerate(a,d,N,p0,20)[1], K)
   plt.errorbar(np.arange(1,N+1)+0.3, m3, s3,fmt='g',alpha=0.5)
   plt.plot(range(1,N+1),crp_particle_num_tables_enumerate(a,d,N,p0, 20, sample=False)[1],'c',linewidth=2)

   plt.plot(range(1,N+1),[expected_num_tables(a, d,i+1, p0)[1] for i in range(N)],'k',linewidth=2)
   plt.plot(range(1, N+1), [crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'b--',linewidth=2)
   plt.grid()
   plt.title(r"Posterior Expected Number of Tables, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(a,d,p0))
   plt.xlabel("# customers of type $s$")
   plt.ylabel('$E[t_s]$')
   plt.xlim(1,N)

def plot_crp_particle_filter_variance(a, d, N, p0, num_particles_list, K=10):
   m1, s1 = get_errorbars(lambda:  [crp_particle_num_tables(a,d,N,p0,i)[1][-1] for i in num_particles_list],K)
   plt.errorbar(num_particles_list, m1, s1,fmt='b',alpha=0.5)
   m2, s2 = get_errorbars(lambda:  [crp_particle_num_tables_enumerate(a,d,N,p0,i)[1][-1] for i in num_particles_list], K)
   plt.errorbar(num_particles_list, m2, s2,fmt='r',alpha=0.5)
   plt.hlines(expected_num_tables(a,d,N,p0)[1],num_particles_list[0], num_particles_list[-1])
   plt.grid()
   plt.title(r"SMC Estimate of $E[t_s]$, $c_s = %d$, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(N,a,d,p0))
   plt.xlabel("# particles")
   plt.ylabel('$E[t_s]$')

def plot_enumerate_particle_filters(a, d, N, p0, num_particles=100):
   #for i in range(20): plt.plot(crp_particle_num_tables_enumerate(a,d,N,p0,100)[1],'g',alpha=0.5)
   before = time.clock()
   plt.plot(crp_particle_num_tables_enumerate(a,d,N,p0,num_particles,sample=False)[1],'b',alpha=0.5)
   plt.plot(crp_particle_num_tables_enumerate(a,d,N,p0,num_particles,sample=True)[1],'r',alpha=0.5)
   plt.plot([crp_fractional_num_tables(a,d,c,p0) for c in range(1,N+1)],'g--')
   after = time.clock()
   print "Elapsed time", after - before
   #for i in range(20): plt.plot(crp_particle_num_tables_enumerate(a,d,N,p0,100, merge=False, sample=False)[1],'r',alpha=0.5)

   before = time.clock()
   plt.plot([expected_num_tables(a, d,i+1, p0)[1] for i in range(N)],'k',linewidth=2)
   after = time.clock()
   print "Elapsed time", after - before
   plt.grid()
   plt.title(r"Posterior Expected Number of Tables, $\alpha=%.1f$, $d=%.1f$, $H(s)=%.1f$"%(a,d,p0))
   plt.xlabel("# customers of type $s$")
   plt.ylabel('# tables')




def main():
  """Create various plots and save them."""
  fig, ax = plot_posterior_components(1., 100)
  fig.set_size_inches(10, 5)
  save_figure(fig, "plot.pdf")

if __name__ == "__main__":
    pass
    #main()
