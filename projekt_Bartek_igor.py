import numpy as np
from wzory_gw import *
f = radians(52)
l = radians(21)
h=100
x,y,z = flh2xyz(f, l, h, a, e2)
print(x,y,z)