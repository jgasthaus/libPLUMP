import libplump
import antaresia.io as aio
import sys
import numpy as np
import scipy.optimize as opt

NUM_FUNC_EVAL = 20

initial = np.array([0, 0.05, 0.7, 0.8, 0.82, 0.84, 0.88, 0.91, 0.92, 0.93, 0.94, 0.95])
#initial = np.array([0.5]*12)

#restaurant = libplump.SimpleFullRestaurant()
#restaurant = libplump.FractionalRestaurant()
restaurant = libplump.HistogramRestaurant()
#restaurant = libplump.KneserNeyRestaurant()
#restaurant = libplump.ReinstantiatingCompactRestaurant()
#restaurant = libplump.StirlingCompactRestaurant()


nodeManager = libplump.SimpleNodeManager(restaurant.getFactory())
parameters = libplump.SimpleParameters()

parameters.discounts = libplump.VectorDouble(initial[1:])
parameters.alpha = initial[0]

train = aio.readBText(sys.argv[1])
valid = aio.readBText(sys.argv[2])
test  = aio.readBText(sys.argv[3])
num_samples = int(sys.argv[4])

#seq = libplump.vectori(range(10))
alldata = train.tolist() + valid.tolist() +test.tolist()
print type(alldata)
print type(alldata[0])
seq = libplump.VectorInt(alldata)
#seq = libplump.VectorInt(map(ord,'oacac'))
#numTypes = max(seq)
numTypes = int(np.max(train) + 1)

model = libplump.HPYPModel(seq, nodeManager, restaurant, parameters, numTypes)
trainLoss = model.computeLosses(0,len(train))

def validLoss(params):
  parameters.discounts = libplump.VectorDouble(params[1:])
  parameters.alpha = params[0]
  probs = model.predictSequence(len(train), len(train)+len(valid))
  loss = np.mean(-np.log2(probs))
  print "evaluated at", params, "loss is", loss
  #print "test loss: ", testLoss()
  return loss
  
def testLoss():
  probs = model.predictSequence(len(train) + len(valid), len(train)+len(valid) + len(test))
  return np.mean(-np.log2(probs))

def testProbs():
  return model.predictSequence(len(train) + len(valid), len(train)+len(valid) + len(test))

def onlineTestLoss():
  losses = model.computeLosses(len(train) + len(valid), len(train)+len(valid) + len(test))
  return np.mean(losses)

bounds = [(0, None)] +  [(0.001, 0.999)]*11

x = initial
probs = []
for i in range(10):
  print "starting optimization iteration", i
  (x, f, d) = opt.fmin_l_bfgs_b(validLoss, x, maxfun=NUM_FUNC_EVAL, bounds = bounds, approx_grad = True)
  print "valid loss", validLoss(x)
  print "test loss ", testLoss()
  probs.append(testProbs())
  print "mean test loss (across samples)", np.mean(-np.log2(np.mean(np.array(probs),0)))
  print "mean test loss (across samples)", np.mean(-np.log2(np.cumsum(np.array(probs)[::-1],0)/(np.arange(len(probs)) + 1)[:,None]))
  for j in range(num_samples):
    print "running sampler"
    model.runGibbsSampler()
  print "valid loss", validLoss(x)
  print "test loss ", testLoss()
print "online test loss", onlineTestLoss()
print "mean loss", np.mean(-np.log2(np.mean(np.array(probs),0)))




# delete model
del model
del nodeManager
del restaurant
