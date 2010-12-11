#!/usr/bin/python
"""Script to demonstrate text prediction with the sequence memoizer."""
import libplump
import sys, tty, termios
import numpy as np

DISCOUNTS = [.62, .69, .74, .80, .95]
CONCENTRATION = 50

def getch():
    """Get a single character from stdin."""
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch


def buildModel(fn):
    """Build a byte-level SM model from the given file."""
    global seq, model, nodeManager, parameters, restaurant
    #restaurant = libplump.SimpleFullRestaurant()
    #restaurant = libplump.HistogramRestaurant()
    #restaurant = libplump.KneserNeyRestaurant()
    restaurant = libplump.ReinstantiatingCompactRestaurant()
    #restaurant = libplump.StirlingCompactRestaurant()
    
    nodeManager = libplump.SimpleNodeManager(restaurant.getFactory())
    parameters = libplump.SimpleParameters(DISCOUNTS, CONCENTRATION)
    
    seq = libplump.VectorInt()
    libplump.pushCharFileToVec(fn, seq)
    numTypes = 256
    
    model = libplump.HPYPModel(seq, nodeManager, restaurant, parameters, numTypes)
    model.computeLosses(0, seq.size())
    
    #model.runGibbsSampler()


def sampleFromModel(prefix='', length=10):
    """Sample symbols from the model."""
    startPos = seq.size()
    for c in prefix:
        seq.push_back(ord(c))
    endPos = seq.size()
    out = []
    for i in xrange(length):
        predictive = model.predictiveDistribution(startPos, endPos)
        sample = libplump.sample_unnormalized_pdf(predictive)
        out.append(chr(sample))
        seq.push_back(sample)
        endPos = seq.size()
    return out


def completion():
    """Simple command line word completion example."""
    key = None
    startPos = seq.size()
    endPos = seq.size()
    while key != 27:
        key = ord(getch())
        if key == 127:
            # backspace
            seq.pop_back()
            endPos -= 1
            sys.stdout.write('\b')
            sys.stdout.write(' ')
            sys.stdout.write('\b')
        elif key == 13:
            #return
            sys.stdout.write('\n')
            startPos = seq.size()
            endPos = seq.size()
        elif key == 9:
            # tab pressed -> predict (map)
            sample = 0
            while chr(int(sample)) != ' ':
                predictive = model.predictiveDistribution(startPos, endPos)
                sample = int(np.argmax(predictive))
                sys.stdout.write(chr(sample))
                seq.push_back(sample)
                endPos += 1
        elif key == 49:
            # 1 pressed -> predict (sample)
            sample = 0
            while chr(int(sample)) != ' ':
                predictive = model.predictiveDistribution(startPos, endPos)
                sample = libplump.sample_unnormalized_pdf(predictive)
                sys.stdout.write(chr(sample))
                seq.push_back(sample)
                endPos += 1
        else:
            sys.stdout.write(chr(key))
            seq.push_back(key)
            endPos += 1
        

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage %s [training_file]" % (sys.argv[0],)
        exit(1)
    print "Building model using particle filtering."
    buildModel(sys.argv[1])
    print "Type text; press <Tab> for MAP completion, 1 for sampling completion,"
    print "<Return> to start a new sequence, <Esc> to exit"
    completion()

    # make sure model is deleted _before_ the underlying nodeManager
    del model
    del nodeManager
