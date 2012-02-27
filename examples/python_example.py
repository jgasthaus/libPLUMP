import libplump

restaurant = libplump.SimpleFullRestaurant()
#restaurant = libplump.HistogramRestaurant()
#restaurant = libplump.KneserNeyRestaurant()
#restaurant = libplump.ReinstantiatingCompactRestaurant()
#restaurant = libplump.StirlingCompactRestaurant()


nodeManager = libplump.SimpleNodeManager(restaurant.getFactory())
parameters = libplump.SimpleParameters()

#seq = libplump.vectori(range(10))
#seq = libplump.vectori([0]*10)
seq = libplump.VectorInt(map(ord,'oacac'))
numTypes = max(seq)

model = libplump.HPYPModel(seq, nodeManager, restaurant, parameters, numTypes)
print model.computeLosses(0,len(seq))
for i in range(seq.size()):
  print model.toString()
  model.runGibbsSampler()

print "Predictions after training:"
for i in range(len(seq)):
   print model.predict(0,i,i)

# save model to file
serializer = libplump.Serializer("model.dump")
serializer.saveNodesAndPayloads(nodeManager, restaurant.getFactory())

# delete model
del model
del nodeManager
del restaurant
del serializer


# load model again from file (note that seq has to be stores separately!)
# set up model
restaurant = libplump.SimpleFullRestaurant()
nodeManager = libplump.SimpleNodeManager(restaurant.getFactory())
model = libplump.HPYPModel(seq, nodeManager, restaurant, parameters, numTypes)
# load structure and parameters from dump
serializer = libplump.Serializer("model.dump")
serializer.loadNodesAndPayloads(nodeManager, restaurant.getFactory())

# this should be the same as before
print "Predictions after save/load:"
for i in range(len(seq)):
   print model.predict(0,i,i)

# make sure destructors are called in correct order (model before nodeManager)
del model
del nodeManager
