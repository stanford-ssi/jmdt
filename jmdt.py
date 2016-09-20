import math
import os

import numpy as np
import matplotlib.pyplot as plt

from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

def run_simulation(# Time step, in seconds.
                   dt = 10,
                   
                   # Send data over to Python over this many integration steps,
                   # i.e. every report_steps*dt seconds.
                   report_steps = 10,
                   
                   # Atmospheric model. 0: None, 1: US1976, 2: NRLMSISE-00
                   # 0 is fastest, 1 is pretty fast (interpolated values),
                   # 2 is much slower but more realistic (changes with time
                   # and space weather).
                   atmosphere = 0,
                   
                   # Earth gravity model. 0: Point mass, 1: WGS84, 2: EGM96
                   # Point mass is fastest, WGS84 has zonal coefficients up to
                   # order 20, EGM96 is painfully slow but includes terms up to
                   # degree 360.
                   earth = 1,
                   
                   # Initial state vector: (x, y, z, vx, vy, vz).
                   state = [6871009, 0, 0, 0, 6620, 3822],
                   
                   # Initial Julian date (sorry).
                   t0 = 2457651.5,

                   # Final time (difference), in seconds.
                   tf = 86400,
                   ):
    process = Popen(['./jmdt'], bufsize=-1,
                        stdout=PIPE, stderr=PIPE, stdin=PIPE)
    
    output_size = 7
    N = math.ceil(tf/dt/report_steps)
    out = np.memmap('output.mmap', dtype=np.double,
                        mode='w+', shape=output_size*N)
    inp = [dt, report_steps, atmosphere, earth]
    inp.extend(state)
    inp.extend([t0, tf, output_size])
    #print '\n'.join(map(str, inp))
    a,b=process.communicate(input='\n'.join(map(str, inp)))
    out.shape = (N, output_size)
    
    ts = out[:, 0]
    xs = out[:, 1]
    plt.plot(ts, xs)
    plt.show()
    print ts
    print repr(a),repr(b)
    
    os.remove('output.mmap')


run_simulation()
