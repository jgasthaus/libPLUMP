"""Script that produces one cell in the log-loss results table (Table 1) in
   the 'Improvements to the Sequence Memoizer' paper.."""

from optparse import OptionParser
import sys
import numpy as np
import libplump as lp

DISCOUNTS = [.62, .69, .74, .80, .95]
NUM_TYPES = 16383
BURN_IN_SAMPLES = 50
PREDICT_SAMPLES = 50

def run(options):
  restaurant = lp.ReinstantiatingCompactRestaurant()
  nodeManager = lp.SimpleNodeManager(restaurant.getFactory())
  discounts = lp.VectorDouble(DISCOUNTS)
  parameters = lp.SimpleParameters(discounts, options.alpha)
  
  seq = lp.VectorInt()
  lp.pushIntFileToVec(options.train_file, seq)
  print >> sys.stderr, "Train seq length: %i" % (seq.size(),)

  # initialize model
  model = lp.HPYPModel(seq, nodeManager, restaurant, parameters, NUM_TYPES)

  #insert training observations into model using particle filter
  model.computeLosses(0, seq.size())

  # add test observations to underlying sequence
  testOffset = seq.size() 
  lp.pushIntFileToVec(options.test_file, seq)
  print >> sys.stderr, "Test seq length: %i" % (seq.size() - testOffset,)

  if options.prediction == 2:
    predictMode = lp.HPYPModel.BELOW
  elif options.prediction == 1:
    predictMode = lp.HPYPModel.FRAGMENT
  else:
    predictMode = lp.HPYPModel.ABOVE

  if options.inference == 1:
    for i in xrange(BURN_IN_SAMPLES):
      print >> sys.stderr, "Burn in iteration %i" %  (i,)
      model.runGibbsSampler()

  if options.prediction != 3:
    loss = float(lp.prob2loss(model.predictSequence(testOffset, seq.size(), predictMode)))
  else:
    loss = float(np.mean(model.computeLosses(testOffset, seq.size())))

  if options.inference == 2 and options.prediction != 3:
    losses = np.zeros((PREDICT_SAMPLES, seq.size() - testOffset))
    for i in xrange(BURN_IN_SAMPLES):
      print >> sys.stderr, "Burn in iteration %i" % (i,)
      model.runGibbsSampler()
    for i in xrange(PREDICT_SAMPLES):
      print >> sys.stderr, "Prediction iteration %i" % (i,)
      model.runGibbsSampler()
      losses[i,:] = model.predictSequence(testOffset, seq.size(), predictMode)
    loss = float(np.mean(-np.log2(np.mean(losses,0))))

  print loss
  
  # make sure destructors are called in correct order
  del model
  del nodeManager

def main():
  parser = OptionParser()
  parser.add_option("--train-file", 
                    dest = "train_file",
                    type = "string",
                    help = "File used for training")
  parser.add_option("--test-file", 
                    dest = "test_file",
                    type = "string", 
                    help = "File used for testing")
  parser.add_option("--alpha",
                    dest = "alpha",
                    type = "float",
                    help = "Concentration parameter",
                    default = 0.0)
  parser.add_option("--inference",
                    dest = "inference",
                    type = "int",
                    help = "Inference mode: 0: particle filter, " + \
                                           "1: 1-sample Gibbs, " + \
                                           "2: 50-sample Gibbs")
  parser.add_option("--prediction",
                    dest = "prediction",
                    type = "int",
                    help = "prediction mode: 0: above, " + \
                                            "1: fragment, " + \
                                            "2: below" + \
                                            "3: particle filter")
  
  (options, args) = parser.parse_args()
  if options.train_file == None or options.test_file == None:
    parser.print_help()
    exit(1)
  if options.inference == 2 and options.prediction == 3:
    print "ERROR: Can't combine particle filter prediction with multiple samples!"
    exit(1)
  run(options)

if __name__ == "__main__":
  main()
