import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE

N = 30

theta = np.linspace(0, np.pi, N)
phi = np.linspace(0, 2*np.pi, N)
print theta
data = np.transpose([np.tile(theta, N), np.repeat(phi, N)])
chunk_size = 200

data = np.array_split(data, np.ceil(float(len(data))/chunk_size))

import os

script = os.path.realpath(__file__).replace('model.py', 'internal.py')

z = []

for x in data:
	process = Popen(['python', script, sys.argv[1], sys.argv[2]], bufsize=-1,
                        stdout=PIPE, stderr=PIPE, stdin=PIPE)
	inp = ""
	for a,b in x:
		inp += "%f,%f " % (a,b)
	a, b = process.communicate(input=inp)
	out = map(float,a.split('\n')[:-1])
	z.extend(out)
	print "extending", len(out)

#z[0] = 42
#z[2] = 42

satname = sys.argv[1].split('.')[0]

#ext = '.solar' if sys.argv[2] == 'solar' else '.drag'
#np.array(z).tofile(satname+ext)
#print "saving", satname+ext
#exit()
z = np.array(z)
z.shape = (N, N)

#import pprint
#print data
#print out

fig=plt.figure()
ax=fig.add_subplot(111)
ax.imshow(z, extent=(0, np.pi, 2*np.pi,0), interpolation='none')
numrows, numcols = z.shape

def format_coord(x, y):
    col = int(x+0.5)
    row = int(y+0.5)
    if col>=0 and col<numcols and row>=0 and row<numrows:
        zz = z[row,col]
        return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, zz)
    else:
        return 'x=%1.4f, y=%1.4f'%(x, y)

ax.format_coord = format_coord
plt.xlabel('theta')
plt.ylabel('phi')
plt.show()

print data
print z

