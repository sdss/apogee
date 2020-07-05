import matplotlib.colors as colors
import matplotlib.cm as cm
import pdb
import numpy as np

def newcmap(low,high) :

    dict = { 'red':   [(0.,0.,0.), (low,0.,0.), (high,1.,1.), (1.,1.,1.)],
             'green': [(0.,0.,0.), (low,0.,0.), (high,1.,1.), (1.,1.,1.)],
             'blue':  [(0.,0.,0.), (low,0.,0.), (high,1.,1.), (1.,1.,1.)] }

    return colors.LinearSegmentedColormap('junk',dict)

def remap(cmap,low,high) :
#    dict = cm.datad[cmap]
    newdict = {}
    if low < 0. : low = 0.
    if low > 1. : low = 1.
    if high < 0. : high = 0.
    if high > 1. : high = 1.
    for color in ['blue','green','red'] :
         new=[]
         new.append((0.,0.,0.))
         new.append((low,0.,0.))
         new.append((high,1.,1.))
         new.append((1.,1.,1.))
         newdict[color] = new

# Python 2
#        b = list(dict[color])
#        new = [b[0]]
#        n=len(b)
#        for i in range(n) :
#            lst = list(b[i])
#            lst[0] = low + lst[0]*(high-low)
#            new.append(tuple(lst))
#        new.append(tuple([1.0,b[n-1][1],b[n-1][2]]))
#         newdict[unicode(color,"utf-8")] = tuple(new)
    return colors.LinearSegmentedColormap('junk',newdict)

