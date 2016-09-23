import math
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from subprocess import Popen, PIPE

import atexit
atexit.register(lambda: os.remove('output.mmap'))

def run_simulation(# Time step, in seconds.
                   dt = 10,
                   
                   # Send data over to Python over this many integration steps,
                   # i.e. every report_steps*dt seconds.
                   report_steps = 1,
                   
                   # Atmospheric model. 0: None, 1: US1976, 2: NRLMSISE-00
                   # 0 is fastest, 1 is pretty fast (interpolated values),
                   # 2 is much slower but more realistic (changes with time
                   # and space weather).
                   atmosphere = 1,
                   
                   # Earth gravity model. 0: Point mass, 1: WGS84, 2: EGM96
                   # Point mass is fastest, WGS84 has zonal coefficients up to
                   # order 20, EGM96 is painfully slow but includes terms up to
                   # degree 360.
                   earth = 1,
                   
                   # Initial state vector: (x, y, z, vx, vy, vz).
                   state = [6871009, 0, 0, 0, 6620, 3822],
                   
                   # Initial Julian date (sorry).
                   t0 = 2457467.50,

                   # Final time (difference), in seconds.
                   tf = 1*86400,
                   
                   # Coefficient of drag.
                   Cd = 2,
                   
                   # Area. Doesn't do anything if satellite loaded from file.
                   A = 0.1*(0.15+0.1+0.15),
                   
                   # Mass of the satellite.
                   mass = 1.0*1.5,
                   
                   # Power simulation parameters. 0: No power simulation,
                   # 1: simulate power (solar panels).
                   power = 1,

                   # Solar panel model file, or "none" if not using one.
                   solar = "ssisat-1/romeo.solar",

                   # Satellite area model file, or "none" if not using one.
                   drag = "ssisat-1/romeo.drag",
                   ):
    #initial = time.time()
    process = Popen(['./jmdt'], bufsize=-1,
                        stdout=PIPE, stderr=PIPE, stdin=PIPE)
    
    output_size = 9
    N = math.floor(tf/dt)/report_steps
    out = np.memmap('output.mmap', dtype=np.double,
                        mode='w+', shape=output_size*N)
    inp = [dt, report_steps, atmosphere, earth]
    inp.extend(state)
    inp.extend([t0, tf, output_size, Cd, A, mass, power, solar, drag])
    print '\n'.join(map(str, inp))
    a,b=process.communicate(input='\n'.join(map(str, inp)))
    print a,b
    
    out.shape = (N, output_size)

    end = np.zeros(shape=(1,output_size))
    while (out[-1] == end).all():
        out = out[:-1]
    
    return out

out = run_simulation()

ts = out[:, 0]
xs = out[:, 1]
ys = out[:, 2]
zs = out[:, 3]
power = out[:, 7]
BC = out[:, 8]

plt.plot(ts, BC)
plt.ylim([0,0.05])
#plt.plot(ts, np.sqrt(xs*xs+ys*ys+zs*zs)/1000.0-6371.009)
#plt.plot(ts, 0*ts)
plt.show()
