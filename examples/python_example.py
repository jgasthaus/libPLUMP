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

for i in range(len(seq)):
   print model.predict(0,i,i)

# make sure destructors are called in correct order
del model
del nodeManager
