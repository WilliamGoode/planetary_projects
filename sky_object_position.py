import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
from orbital import planetPlot

print('Sky object positions development test...')

# Load spice kernels:
spice.furnsh('metakernel_pc.txt')

# Astro parameters:
r_earth = spice.bodvrd('EARTH', 'RADII', 3)[1] # [km]
re = r_earth[0]
rp = r_earth[2]
f = (re -rp)/re
print('Earth radii:', r_earth)
AU = 1.4959787e8 # [km]
day = 86400 # [s]
moon_sma = 384399 # [km]

date_utc = '2024 March 2, 14:00:00 UTC' 
date_et = spice.str2et(date_utc)
print('ET: ', date_et)

et_array = np.linspace(date_et, date_et+27*day, 500)

moon_states = spice.spkpos('MOON', et_array, 'ECLIPJ2000', 'none', 'EARTH')
# print(moon_states)
moon_pos = moon_states[0] # intertial
# moon_pos = moon_pos.T
# print(moon_pos)

# Convert positions to ECEF frame:
# ECI2ECEF = spice.pxform('ECLIPJ2000', 'IAU_EARTH', date_et)
# print('ECI to ECEF matrix:')
# print(ECI2ECEF)
# moon_pos_ecef = spice.mxvg(ECI2ECEF, moon_pos)
# print('Moon pos in ECEF:')
# print(moon_pos_ecef)

# plt.figure()
# plt.plot(moon_pos[:,0], moon_pos[:,1], label='Moon Position')
# plt.xlabel('x (km)')
# plt.ylabel('y (km)')
# plt.xlim([-moon_sma, moon_sma])
# plt.ylim([-moon_sma, moon_sma])
# plt.legend()
# plt.grid()
# plt.show()

X_earth, Y_earth, Z_earth = planetPlot(re)

# Plotting Earth and Orbit
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X_earth, Y_earth, Z_earth, color='blue', alpha=0.7)
ax.plot3D(moon_pos[:,0], moon_pos[:,1], moon_pos[:,2], 'black')
ax.scatter(moon_pos[-1,0], moon_pos[-1,1], moon_pos[-1,2], c='black')
ax.view_init(90, 0)  # Changing viewing angle (adjust as needed)
plt.title('Lunar Position')
ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')
# Make axes limits
# xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()]).T
# xyzlim = np.array([5*r_earth, 5*r_earth, 5*r_earth]).T
# XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
# ax.set_xlim3d(XYZlim)
# ax.set_ylim3d(XYZlim)
# ax.set_zlim3d(XYZlim)
limit = 1.0*moon_sma
ax.axes.set_xlim3d(left=-limit, right=limit) 
ax.axes.set_ylim3d(bottom=-limit, top=limit) 
ax.axes.set_zlim3d(bottom=-limit, top=limit)
plt.show()
