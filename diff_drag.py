import math
import os
import time

import matplotlib as mpl
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.basemap import Basemap

from subprocess import Popen, PIPE

import atexit
#atexit.register(lambda: os.remove('output.mmap'))

# Radius of the Earth
r_earth = 6371009 # [m]

# Initial altitude
altitude_0 = 600 * 1000 # [m]

# TODO: Currently assumes orbit starts at periapsis
# Periapsis
r_peri = altitude_0 + r_earth

# TODO: Currently assumes orbit starts at periapsis
# Orbital elements
orbital_e = 0.0                                                 # Eccentricity [unitless]
orbital_i = 28.0                                                    # Inclination [deg]
orbital_a = r_peri * (1 + ((1 + orbital_e)/(1 - orbital_e)))/2.0    # Semi-major axis [m]
#orbital_Omega                                                      # unimplemented
#orbital_omega                                                      # unimplemented
#orbital_nu                                                         # unimplemented

# Gravitational constant
mu_earth = 398600441500000.0 # [m^3/s^2] - definitely not that exact to that many decimals, but converted from km to m

# Calculate v_0 from orbital elements (see Vis-Viva Law)
v_0 = (mu_earth * ((2.0/r_peri)-(1.0/orbital_a))) ** 0.5    # [m/s]

# TODO: Currently assumes orbit starts at periapsis and that periapsis occurs when x = 0
# Define three components of initial velocity
vx_0 = 0.0                              # [m/s]
vy_0 = v_0 * math.sin(orbital_i)        # [m/s]
vz_0 = v_0 * math.cos(orbital_i)        # [m/s]

# Calculate initial velocity magnitude
# v_0 = np.sqrt(vx_0**2 + vy_0**2 + vz_0**2)

# Specify desired separation velocity between satellites (set to 0 if not separating or using a single satellite)
vseparation = 0.1 # [m/s]

# Specify desired separation distance
target_distance = 100000;

vx_0_1 = vx_0 + ((vx_0/v_0) * vseparation/2)
vy_0_1 = vy_0 + ((vy_0/v_0) * vseparation/2)
vz_0_1 = vz_0 + ((vz_0/v_0) * vseparation/2)

vx_0_2 = vx_0 - ((vx_0/v_0) * vseparation/2)
vy_0_2 = vy_0 - ((vy_0/v_0) * vseparation/2)
vz_0_2 = vz_0 - ((vz_0/v_0) * vseparation/2)

def run_simulation(# Time step, in seconds.
                   dt = 10,

                   # Send data over to Python over this many integration steps,
                   # i.e. every report_steps*dt seconds.
                   report_steps = 10,

                   # Atmospheric model. 0: None, 1: US1976, 2: NRLMSISE-00
                   # 0 is fastest, 1 is pretty fast (interpolated values),
                   # 2 is much slower but more realistic (changes with time
                   # and space weather).
                   atmosphere = 2,

                   # Earth gravity model. 0: Point mass, 1: WGS84, 2: EGM96
                   # Point mass is fastest, WGS84 has zonal coefficients up to
                   # order 20, EGM96 is painfully slow but includes terms up to
                   # degree 360.
                   earth = 1,

                   # TODO: Currently assumes orbit starts at periapsis and that periapsis occurs when x = 0
                   # Initial state vector: (x, y, z, vx, vy, vz).
                   state = [r_peri, 0, 0, vx_0_1, vy_0_1, vz_0_1],

                   # Initial Julian date (sorry).
                   t0 = 2457467.50, # [days] - yes, it's dumb, but it's standard

                   # Final time (difference)
                   tf = 30 * 86400, # [s]

                   # Coefficient of drag.
                   Cd = 2,

                   # Area. Doesn't do anything if satellite loaded from file.
                   A = 0.1*(0.15+0.1+0.15),

                   # Mass of the satellite.
                   mass = 4.0, # [kg]

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

                   # Separation target
                   separation_target = target_distance, # [m]

                   # Initial state of the second satellite.
                   second_state = [r_earth + altitude_0, 0, 0, vx_0_2, vy_0_2, vz_0_2],

                   # Orientation mode for the first satellite.
                   # t: Target (i.e. the other satellite)
                   # z: Zenith
                   # n: Nadir
                   # p: Prograde
                   # r: Retrograde
                   # s: Sun-facing
                   first_orientation = "r",

                   # Orientation mode for the second satellite.
                   second_orientation = "z",

                   # Constants of PI controller
                   k_i = -0.02,
                   k_p = -7500.0

                   ):
    #initial = time.time()
    process = Popen(['./jmdt'], bufsize=-1,
                        stdout=PIPE, stderr=PIPE, stdin=PIPE)

    output_size = 21
    N = math.floor(tf/dt)/report_steps
    out = np.memmap('output.mmap', dtype=np.double,
                        mode='w+', shape=output_size*N)
    inp = [dt, report_steps, atmosphere, earth]
    inp.extend(state)
    inp.extend([t0, tf, output_size, Cd, A, mass, power, solar, drag,
                solar_efficiency, two_satellites, separation_target])
    inp.extend(second_state)
    inp.extend([first_orientation, second_orientation])
    inp.extend([k_i, k_p])
    #print '\n'.join(map(str, inp))
    a,b=process.communicate(input='\n'.join(map(str, inp)))
    print a,b

    out.shape = (N, output_size)

    end = np.zeros(shape=(1,output_size))
    while (out[-1] == end).all():
        out = out[:-1]

    return out

out = run_simulation()

ts = out[:, 0]/86400.0
xs = out[:, 1]
ys = out[:, 2]
zs = out[:, 3]
xv = out[:, 4]
yv = out[:, 5]
zv = out[:, 6]
power = out[:, 7]
BC = out[:, 8]
xs2 = out[:, 9]
ys2 = out[:, 10]
zs2 = out[:, 11]
xv2 = out[:, 12]
yv2 = out[:, 13]
zv2 = out[:, 14]

dd = out[:, 17]
dd2 = out[:, 18]

#plt.plot(ts, power)

fig1 = plt.figure(figsize = (8,6))

# # Separation distance plot
# ax_pos = plt.subplot(311)
#
# plt.plot(ts, np.sqrt((xs-xs2)**2 + (ys-ys2)**2 + (zs-zs2)**2)/1000.0 , 'g-')
# plt.ylabel('Separation distance (km)')
# plt.xlabel('Time (days)')
# plt.savefig('diff_drag.png', bbox_inches='tight')

# Altitude plot
ax_pos = plt.subplot(211)
plt.title('Differential Drag Performance in 600 x 600 km 28-degree Orbit')
plt.plot(ts, np.sqrt(xs*xs + ys*ys + zs*zs)/1000.0 - 6371.009, 'r-', ts, np.sqrt(xs*xs + ys*ys + zs*zs)/1000.0 - 6371.009, 'b-')
plt.ylabel('Satellite altitude (km)')

# Velocity plot
ax_vel = plt.subplot(212)
plt.plot(ts, np.sqrt(xv*xv + yv*yv + zv*zv))
plt.ylabel('Satellite velocity (m/s)')
plt.xlabel('Time (days)')
plt.savefig('diff_drag.png', bbox_inches='tight')

fig2 = plt.figure(figsize = (8,6))
ax_dis = plt.subplot(211)
plt.title('Differential Drag Performance in 600 x 600 km 28-degree Orbit')
plt.plot(ts, np.sqrt((xs-xs2)**2 + (ys-ys2)**2 + (zs-zs2)**2)/1000.0 , 'g-')
plt.plot(ts, np.ones(len(ts)) * target_distance/1000 , 'k--')
plt.ylabel('Separation Distance (km)')

# ax_dd = plt.subplot(312)
# plt.plot(ts, dd)
# plt.ylabel('Drift Velocity (m/s)')

ax_dd2 = plt.subplot(212)
plt.plot(ts, dd2)
plt.ylabel('Filtered Drift Velocity (m/s)')
plt.xlabel('Time (days)')
plt.savefig('pi_controller.png', bbox_inches='tight')

plt.show()

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(xs, ys, zs, color='#be1e2d')

#
# # set up orthographic map projection with
# # perspective of satellite looking down at 50N, 100W.
# # use low resolution coastlines.
# map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# # draw coastlines, country boundaries, fill continents.
# map.drawcoastlines(linewidth=0.25)
# map.drawcountries(linewidth=0.25)
# map.fillcontinents(color='coral',lake_color='aqua')
# # draw the edge of the map projection region (the projection limb)
# map.drawmapboundary(fill_color='aqua')
# # draw lat/lon grid lines every 30 degrees.
# map.drawmeridians(np.arange(0,360,30))
# map.drawparallels(np.arange(-90,90,30))
# # make up some data on a regular lat/lon grid.
# nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
# lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
# lons = (delta*np.indices((nlats,nlons))[1,:,:])
# wave = 0.75*(np.sin(2.*lats)**8*np.cos(4.*lons))
# mean = 0.5*np.cos(2.*lats)*((np.sin(2.*lats))**2 + 2.)
# # compute native map projection coordinates of lat/lon grid.
# x, y = map(lons*180./np.pi, lats*180./np.pi)
# # contour data over the map.
# cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
# plt.title('contour lines over filled continent background')
# plt.show()



#plt.ylim([0,0.05])
#plt.plot(ts, np.sqrt(xs*xs+ys*ys+zs*zs)/1000.0-6371.009)
#plt.plot(ts, 0*ts)
