import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
from orbital import planetPlot
import orbital

print('Sky object positions development test...')

# Load spice kernels:
mk_file = '../SPICE_kernels/mk/metakernel_mac.txt'
spice.furnsh(mk_file)

# Astro parameters:
r_earth = spice.bodvrd('EARTH', 'RADII', 3)[1] # [km]
re = r_earth[0]
rp = r_earth[2]
f = (re -rp)/re
print('Earth radii:', r_earth)
r_moon = spice.bodvrd('MOON', 'RADII', 3)[1] # [km]
re_moon = r_moon[0]
AU = 1.4959787e8 # [km]
day = 86400 # [s]
moon_sma = 384399 # [km]

date_utc = '2024 March 17, 21:03:00 UTC' 
date_et = spice.str2et(date_utc)
print('ET: ', date_et)

et_array = np.linspace(date_et, date_et+1*day, 100)

moon_states = spice.spkpos('MOON', et_array, 'ECLIPJ2000', 'none', 'EARTH')
# print(moon_states)
moon_pos = moon_states[0] # intertial
# moon_pos = moon_pos.T

# Convert positions to ECEF frame:
moon_pos_ecef = np.empty((moon_pos.shape[0],moon_pos.shape[1]))
for i in range(moon_pos.shape[0]):
    ECI2ECEF = spice.pxform('ECLIPJ2000', 'IAU_EARTH', et_array[i])
    moon_pos_ecef[i,:] = spice.mxvg(ECI2ECEF, moon_pos[i,:])

# Get observation location in ECEF frame:
lat_test = 39.5*(np.pi/180) # [rad]
lon_test = -104.7*(np.pi/180) # [rad]
alt_test = 1828e-3 # [km]

moon_pos_ecef_test = moon_pos_ecef[1,:]

obs_pos_ecef = spice.georec(lon_test, lat_test, alt_test, re, f)
los_ecef = moon_pos_ecef_test - obs_pos_ecef

pos_enu = orbital.C_ecef2enu(lat_test, lon_test)@los_ecef
print('Moon position in ENU:',pos_enu)
az, el = orbital.AZEL(pos_enu)
print('Az (deg):', az*180/np.pi)
print('El (deg):', el*180/np.pi)


fig = plt.figure()
ax = fig.add_subplot(polar=True)
ax.scatter(az, el*180/np.pi, c='red', s=50)
ax.set_theta_zero_location("N")  # theta=0 at the top
ax.set_theta_direction(-1)  # theta increasing clockwise
ax.set_rlim(bottom=90, top=0)
# plt.show()

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
X_moon, Y_moon, Z_moon = planetPlot(re_moon)

# Plotting Earth and Orbit
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X_earth, Y_earth, Z_earth, color='blue', alpha=0.7)
ax.plot3D(moon_pos[:,0], moon_pos[:,1], moon_pos[:,2], 'black', linewidth=0.5)
ax.scatter(moon_pos[-1,0], moon_pos[-1,1], moon_pos[-1,2], c='red', s=2)
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
