import math
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

from subprocess import Popen, PIPE

import atexit
#atexit.register(lambda: os.remove('output.mmap'))

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
                   tf = 10*86400,
                   
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
                   solar = "ssisat-1/rev2.solar",

                   # Satellite area model file, or "none" if not using one.
                   drag = "ssisat-1/rev2.drag",
                   
                   # Efficiency of the solar panels.
                   solar_efficiency = 0.27,
                   
                   # Simulate two satellites at the same time. No
                   # multithreading for now, though.
                   two_satellites = 1,
                   
                   # Initial state of the second satellite.
                   #second_state = [6792969.91294023, 896908.98632624, 517817.15727738, -1145.41077485, 6544.81774577, 3778.47253018], #[6871009, 0, 0, 0, 6620, 3822],
                   #second_state = [6870294.57853644, 86057.01722171, 49684.27308189, -109.90904956, 6619.31167374, 3821.6014847],
                   second_state = [6871009, 0, 0, 0, 6620, 3822],
                   
                   # Orientation mode for the first satellite.
                   # t: Target (i.e. the other satellite)
                   # r: Radial
                   # p: Prograde
                   # s: Sun-facing
                   first_orientation = "c1",
                   
                   # Orientation mode for the second satellite.
                   second_orientation = "c2"
                   ):
    #initial = time.time()
    process = Popen(['./jmdt'], bufsize=-1,
                        stdout=PIPE, stderr=PIPE, stdin=PIPE)
    
    output_size = 17
    N = math.floor(tf/dt)/report_steps
    out = np.memmap('output.mmap', dtype=np.double,
                        mode='w+', shape=output_size*N)
    inp = [dt, report_steps, atmosphere, earth]
    inp.extend(state)
    inp.extend([t0, tf, output_size, Cd, A, mass, power, solar, drag,
                solar_efficiency, two_satellites])
    inp.extend(second_state)
    inp.extend([first_orientation, second_orientation])
    print '\n'.join(map(str, inp))
    a,b=process.communicate(input='\n'.join(map(str, inp)))
    print a,b
    
    out.shape = (N, output_size)
    print out

    end = np.zeros(shape=(1,output_size))
    while (out[-1] == end).all():
        out = out[:-1]
    
    return out

out = run_simulation()

ts = out[:, 0]
xs = out[:, 1]
ys = out[:, 2]
zs = out[:, 3]

vx = out[:, 4]
vy = out[:, 5]
vz = out[:, 6]
rdots = np.array([vx, vy, vz]).T

rs = np.array([xs, ys, zs]).T
power = out[:, 7]
BC = out[:, 8]

x2 = out[:, 9]
y2 = out[:, 10]
z2 = out[:, 11] 
r2 = np.array([x2, y2, z2]).T

r0 = rs[0]
BC2 = out[:, 16]

"""
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

R = 6371009
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
x=R*np.cos(u)*np.sin(v)
y=R*np.sin(u)*np.sin(v)
z=R*np.cos(v)
ax.plot_wireframe(x, y, z, color="r")
for i in range(10):
    ax.scatter(xs[i], ys[i], zs[i], color='blue')
for i in range(10):
    ax.scatter(x2[i], y2[i], z2[i], color='black')
plt.show()
exit()"""

#for i in range(len(rs)):
#    print ts[i], np.linalg.norm(rs[i]-r0), rs[i], rdots[i]
#    pass#print np.linalg.norm(r-r0)

day = 86400.
#plt.plot(ts, np.linalg.norm(rs-r0, axis=1))
#plt.plot(ts/day, BC)
#plt.plot(ts/day, BC2)
plt.plot(ts/day, np.linalg.norm(rs-r2, axis=1))
plt.plot(ts/day, ts*0+100000)
#plt.ylim([0, 200000])
#plt.ylim([0,0.05])
#plt.plot(ts, np.sqrt(xs*xs+ys*ys+zs*zs)/1000.0-6371.009)
#plt.plot(ts, 0*ts)
plt.show()
