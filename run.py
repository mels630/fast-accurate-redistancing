from matplotlib import pyplot as pyp
import pyRedist

pyr = pyRedist.pyRedist()
print "Two-dimensional checks"
print "Checking fast-marching implementation"
for pw in xrange(5,9):
    pyr.Run(2, 2**pw, 2**pw, 1)
print "Checking quadratic convergence"
for pw in xrange(5,9):
    pyr.Run(2, 2**pw, 2**pw, 2)
print "Checking cubic convergence"
for pw in xrange(5,9):
    pyr.Run(2, 2**pw, 2**pw, 3)

print "Three-dimensional checks"
print "Checking fast-marching implementation"
for pw in xrange(4,7):
    pyr.Run(3, 2**pw, 2**pw, 1)
print "Checking quadratic convergence"
for pw in xrange(4,7):
    pyr.Run(3, 2**pw, 2**pw, 2)
