
import numpy as np
import matplotlib.pyplot as plt
import spiceypy as sp
import scipy.io

sp.furnsh('metakernel_pc.txt')
orbits = scipy.io.loadmat('planet_orbits.mat')

planet_names = ["Earth", "Mars"]

date_1_utc = '1 Aug 1989'
date_2_utc = '1 Jul 1990'

AU = 1.4959787e8 # [km]
day = 86400 # [s]
a_earth = AU
a_mars = 1.523679342 * AU

frame = 'ECLIPJ2000'

date_1_et = sp.str2et(date_1_utc)
date_2_et = sp.str2et(date_2_utc)

et_array = np.arange(date_1_et, date_2_et, day)

numberOfPlanets = len(planet_names)

earth_pos, lightTime = sp.spkpos('EARTH', et_array, frame, 'NONE', 'SOLAR SYSTEM BARYCENTER')

earth_pos = earth_pos.T

earth_pos_au = np.divide(earth_pos, AU)

fig = plt.Figure()
plt.plot(earth_pos_au[0], earth_pos_au[1])
plt.show()

print('Stop here')