# Orbital Motion python tools
# Author: William Goode
# Date: 28 Jan 2023

# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import spiceypy as spice

class OrbitalElements(object):
    a = None
    e = None
    i = None
    Omega = None
    omega = None
    f = None

# Euler rotation matices
def rot1(alpha):
    return np.array([[1., 0., 0.], 
                     [0., np.cos(alpha), np.sin(alpha)], 
                     [0., -np.sin(alpha), np.cos(alpha)]])

def rot2(alpha):
    return np.array([[np.cos(alpha), 0., -np.sin(alpha)], 
                     [0., 1., 0.], 
                     [np.sin(alpha), 0., np.cos(alpha)]])

def rot3(alpha):
    return np.array([[np.cos(alpha), np.sin(alpha), 0.], 
                     [-np.sin(alpha), np.cos(alpha), 0.], 
                     [0., 0., 1.]])

def C_ecef2enu(lat,lon):
    return np.array([[-np.sin(lon), np.cos(lon), 0],
                     [-np.cos(lon)*np.sin(lat), -np.sin(lon)*np.sin(lat), np.cos(lat)],
                     [np.cos(lon)*np.cos(lat), np.sin(lon)*np.cos(lat), np.sin(lat)]])

def AZEL(pos_enu):
    az = np.arctan2(pos_enu[0],pos_enu[1])
    if az < 0:
        az = az + 2*np.pi
    el = np.arcsin(pos_enu[2]/np.linalg.norm(pos_enu))
    return az, el

def gregorian_to_julian(year, month, day):
    i = int((month - 14) / 12)
    jd = day - 32075
    jd += int((1461 * (year + 4800 + i)) / 4)
    jd += int((367 * (month - 2 - (12 * i))) / 12)
    jd -= int((3 * int((year + 4900 + i) / 100)) / 4)
    return jd

def oscelts2state(elements, mu):
    """
    Get state [x, y, z, xdot, ydot, zdot] from orbital elements.
    """
    # a, e, i, Omega, omega, f
    a = elements.a
    e = elements.e
    i = elements.i
    Omega = elements.Omega
    omega = elements.omega
    f = elements.f
    p = a*(1-e**2)
    if e < 1e-11 and i < 1e-11: # circular equatorial
        omega = 0.0
        Omega = 0.0
    elif e < 1e-11: # circular inclined
        omega = 0.0
        f = 0.0
    elif abs(i) < 1e-11: # elliptical equatorial
        Omega = 0.0
        # omega = omega_true
    
    r = p/(1+e*np.cos(f)) # radial distance [km]
    # Vectors in the perifocal (pqw) coordinate system:
    r_pqw = np.array([r*np.cos(f), r*np.sin(f), 0.0])
    v_pqw = np.array([-np.sqrt(mu/p)*np.sin(f), np.sqrt(mu/p)*(e + np.cos(f)), 0.0])

    # Rotation matrices from pqw to ijk:
    R1 = rot3(-omega)
    R2 = rot1(-i)
    R3 = rot3(-Omega)

    r_ijk = R3@R2@R1@r_pqw
    v_ijk = R3@R2@R1@v_pqw

    state = np.concatenate([r_ijk, v_ijk], axis=None)
    # print(r_ijk)
    # print(v_ijk)
    # print(state)

    return state
    
def state2oscelt(state, mu, et):

    """
    Get orbital elements (a, e, i, Omega, omega, f) from state [x, y, z, xdot, ydot, zdot]
    """

    r_vec = state[0:3]
    v_vec = state[3:6]
    v = np.linalg.norm(v_vec)
    r = np.linalg.norm(r_vec)

    elements = OrbitalElements()

    epsilon = 1e-11 # threshold for special cases

    i_hat = np.array([1.,0.,0.])
    j_hat = np.array([0.,1.,0.])
    k_hat = np.array([0.,0.,1.])

    h_vec = np.cross(r_vec,v_vec)
    h = np.linalg.norm(h_vec)
    n_vec = np.cross(k_hat,h_vec)
    n = np.linalg.norm(n_vec)
    e_vec = ((v**2 - mu/r)*r_vec - np.dot(r_vec,v_vec)*v_vec)/mu
    e = np.linalg.norm(e_vec)
    energy = 0.5*v**2 - mu/r
    if e != 1.0:
        a = -mu/(2*energy)
        p = a*(1-e**2)
    else:
        p = h**2/mu
        a = -99

    elements.a = a
    elements.e = e

    cos_i = h_vec[2]/h

    if e > epsilon:
        cos_f = np.dot(e_vec, r_vec)/(e*r)
        f = np.arccos(cos_f) # [rad]
    else:
        f = 0.0
    
    inc = np.arccos(cos_i) # [rad]
    if abs(inc) > epsilon: # is the orbit inclined?
        cos_LAN = n_vec[0]/n
        cos_AP = np.dot(n_vec, e_vec)/(n*e)
        AP = np.arccos(cos_AP) # [rad]
    else:
        AP = 0.0

    if e > epsilon and abs(inc) > epsilon: # Elliptical inclined
        print('Elliptical inclined')
        LAN = np.arccos(cos_LAN) # [rad]
        if n_vec[1] < 0:
            LAN = 2*np.pi - LAN
        if e_vec[2] < 0:
            AP = 2*np.pi - AP
        if np.dot(r_vec,v_vec) < 0:
            f = 2*np.pi - f
        elements.Omega = LAN
        elements.omega = AP
        elements.f = f
    # Special Cases:
    elif e > epsilon and abs(inc) < epsilon: # Elliptical equatorial
        print('Elliptical equatorial')
        cos_APtrue = e_vec[0]/e
        APtrue = np.arccos(cos_APtrue) # [rad]
        if e_vec[1] < 0:
            APtrue = 2*np.pi - APtrue
        elements.omega = APtrue
    elif e < epsilon and abs(inc) > epsilon: # Circular Inclined
        print('Circular Inclined')
        cos_u = np.dot(n_vec,r_vec)/(n*r)
        u = np.arccos(cos_u) # [rad]
        if r_vec[2] < 0:
            u = 2*np.pi - u
    elif e < epsilon and abs(inc) < epsilon: # Circular Equatorial
        print('Circular Equatorial')
        cos_lamtrue = r_vec[0]/r
        lamtrue = np.arccos(cos_lamtrue) # [rad]
        if r_vec[1] < 0:
            lamtrue = 2*np.pi - lamtrue
        elements.omega = 0.0
        elements.Omega = 0.0
        elements.f = 0.0
    else:
        print('Orbit could not be determined.')

    # print('inc', inc, ' radians')
    # elements.f = f
    elements.i = inc
    # elements.Omega = LAN
    # elements.omega = AP

    return elements


def numintegrate2body(state0, t, mu):

    """
    Propogate orbit on 2-body problem given initial state [x, y, z, xdot, ydot, zdot].
    """

    X0 = state0 # initial conditions
    
    def dXdT(X, t):
        vx = X[3]
        vy = X[4]
        vz = X[5]

        r3 = X[0]**2 + X[1]**2 + X[2]**2
        xdotdot = (-mu/r3**(3/2))*X[0]
        ydotdot = (-mu/r3**(3/2))*X[1]
        zdotdot = (-mu/r3**(3/2))*X[2]

        dX = [vx, vy, vz, xdotdot, ydotdot, zdotdot]
        return dX

    sol = odeint(dXdT, state0, t)
    r_sat = [sol[:,0], sol[:,1], sol[:,2]]
    v_sat = [sol[:,3], sol[:,4], sol[:,5]]

    state = np.transpose(sol)
    return state

def energy(state, mu): # specific energy
    vsquared = (state[3]**2 + state[4]**2 + state[5]**2)
    r = (state[0]**2 + state[1]**2 + state[2]**2)**(1/2)
    return (1/2)*vsquared - mu/r

def planetPlot(r_planet):
    # Setting up Spherical Earth to Plot
    N = 50
    phi = np.linspace(0, 2 * np.pi, N)
    theta = np.linspace(0, np.pi, N)
    theta, phi = np.meshgrid(theta, phi)

    X_planet = r_planet * np.cos(phi) * np.sin(theta)
    Y_planet = r_planet * np.sin(phi) * np.sin(theta)
    Z_planet = r_planet * np.cos(theta)

    return X_planet, Y_planet, Z_planet
    
def main():

    # lsk_file = 'naif0012.tls'
    # spice.furnsh(lsk_file)
    mk_file = 'metakernel_pc.txt'
    spice.furnsh(mk_file)

    utc =  '2004-06-11T19:32:00'
    et = spice.str2et(utc)
    print(et)

    print("Testing orbital functions...")
    r_earth = 6378.14  # Average radius of Earth [km]
    X_earth, Y_earth, Z_earth = planetPlot(r_earth)

    mu_earth = 3.986004418E+05  # Earth's gravitational parameter  [km^3/s^2]

    time_arr = np.linspace(0, 6*3600, 500)

    # # Initial Conditions for eccentric inclined orbit
    # X_0 = -2500  # [km]
    # Y_0 = -5500  # [km]
    # Z_0 = 3400  # [km]
    # VX_0 = 7.5  # [km/s]
    # VY_0 = 0.0  # [km/s]
    # VZ_0 = 4.0  # [km/s]
    # state_0 = np.array([X_0, Y_0, Z_0, VX_0, VY_0, VZ_0])
    # print('Initial state:', state_0)
    # print(state_0.shape)

    # Initial Conditions for circulalar equatorial orbit
    X_0 = r_earth + 850  # [km]
    Y_0 = 0.0  # [km]
    Z_0 = 0.0  # [km]
    VX_0 = 0.0  # [km/s]
    VY_0 = np.sqrt(mu_earth/np.linalg.norm(X_0))  # [km/s]
    VZ_0 = 0.0  # [km/s]
    state_0 = np.array([X_0, Y_0, Z_0, VX_0, VY_0, VZ_0])
    print('Initial state:', state_0)
    print(state_0.shape)

    # Convert state to elements:
    oscelts = state2oscelt(state_0, mu_earth, et)
    print('a, e, i, Omega, w, f:', oscelts.a, oscelts.e, oscelts.i, oscelts.Omega, oscelts.omega, oscelts.f)

    # Test elements to create orbit:

    # # Convert elements to state:
    state2 = oscelts2state(oscelts, mu_earth)
    print('state 2: ', state2)
    print(state2.shape)

    # Modify elemtents:
    oscelts.i = 30.0 * (np.pi/180) # degrees
    state3 = oscelts2state(oscelts, mu_earth)

    state = numintegrate2body(state_0, time_arr, mu_earth)
    state_2 = numintegrate2body(state2, time_arr, mu_earth) # run the new initial state from elemtents
    state_3 = numintegrate2body(state3, time_arr, mu_earth) # run modified elements of orbit
    test_state = state[:,0:4]
    print(test_state.shape)
    
    # oscelts = state2oscelt(test_state, mu_earth, et)

    # Plotting Earth and Orbit
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X_earth, Y_earth, Z_earth, color='blue', alpha=0.7)
    ax.plot3D(state[0,:], state[1,:], state[2,:], 'black')
    ax.plot3D(state_2[0,:], state_2[1,:], state_2[2,:], 'blue')
    ax.plot3D(state_3[0,:], state_3[1,:], state_3[2,:], 'red')
    ax.view_init(90, 0)  # Changing viewing angle (adjust as needed)
    plt.title('Two-Body Orbit')
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
    limit = 1.5*r_earth
    ax.axes.set_xlim3d(left=-limit, right=limit) 
    ax.axes.set_ylim3d(bottom=-limit, top=limit) 
    ax.axes.set_zlim3d(bottom=-limit, top=limit)
    plt.show()

if __name__ == '__main__':
    main()