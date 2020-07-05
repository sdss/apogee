# David Nidever's Python Utility Functions

This package has a bunch of small functions that I find useful while working in python.
Most of the functions are in a module called "dlnpyutils".

# Installation

Until I get this in pypi the best way is the usual git clone and setup install:

```
git clone https://github.com/dnidever/dlnpyutils.git
cd dlnpyutils
python setup.py install
```

# Using the package

To import all of the package functions into the namespace do:
```python
from dlnpyutils.dlnpyutils import *
```

# Some of the functions

 dlnpyutils:
 - mad: median absolute deviation of array
 - minmax: minimum and maximum of an array
 - stat: many useful statistics of an array
 - strlen: number of characters in a string array or list
 - strip: strip whitespace from string array or list
 - strjoin: combine string arrays or scalars
 - strsplit: split string arrays
 - pathjoin: join two pathname components
 - first_el: return the first element of an array or list
 - grep: grep on a string array
 - readlines: read a file into a string array
 - writelines: write a string array to a file
 - remove_indices: remove certain indices from an array
 - numlines: return the number of lines in a file
 - basiclogger: return a basic logger to the screen and optionally a file
 - remove: delete multiple files and allow for non-existence
 - lt: takes the lesser of x or limit
 - gt: takes the greater of x or limit
 - limit: require x to be within upper and lower limits
 - gaussian: return Gaussian plus constant
 - gaussfit: fit a 1-D Gaussian to X/Y data
 - poly: evaluate a polynomial function of a variable
 - slope: derivative or slope of an array

 job_daemon:
 This is a simple batch job manager.
