"""Worker script that produces one cell in the log-loss results table."""

from optparse import OptionParser, OptionGroup
import sys
from numpy import mean
from antaresia.filecollection import FileCollection
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
    loss = lp.prob2loss(model.predictSequence(testOffset, seq.size(), predictMode))
  else:
    loss = mean(model.computeLosses(testOffset, seq.size()))

  if options.inference == 2 and options.prediction != 3:
    losses = []
    for i in xrange(BURN_IN_SAMPLES):
      print >> sys.stderr, "Burn in iteration %i" % (i,)
      model.runGibbsSampler()
    for i in xrange(PREDICT_SAMPLES):
      print >> sys.stderr, "Prediction iteration %i" % (i,)
      model.runGibbsSampler()
      losses.append(lp.prob2loss(model.predictSequence(testOffset, seq.size(), predictMode)))
    loss = mean(losses)

  print loss
  fc = FileCollection(options.collection)
  doc = vars(options)
  doc["loss"] = loss
  fc.insert(doc)
  
  # make sure destructors are called in correct order
  del model
  del nodeManager

def main():
  parser = OptionParser()
  parser.add_option("--train-file", 
                    dest = "train_file",
                    type = "string")
  parser.add_option("--test-file", 
                    dest = "test_file",
                    type = "string")
  parser.add_option("--collection", 
                    dest = "collection",
                    type = "string",
                    default="collections/default")
  parser.add_option("--alpha",
                    dest = "alpha",
                    type = "float")
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
  if options.inference == 2 and options.prediction == 3:
    exit(1)
  run(options)

if __name__ == "__main__":
  main()
