import numpy as np
import matplotlib.pyplot as plt

# Parameters
apogee = 3000.0           # m

def air_density(h):
    # Constants
    t0 = 288.15         # Sea-level temperature in K
    pho0 = 1.225        # Sea-level density in kg/m³
    r = 287.058         # Specific gas constant in J/(kg·K)
    P = 101325          # Sea-level pressure in Pa
    lapse_rate = 0.0065 # Lapse rate in K/m
    g = 9.81            # Acceleration due to gravity in m/s²

    T = t0 - lapse_rate * h
    pressure = P * (T / t0) ** (g / (lapse_rate * r))
    density = pressure / (r * T)
    return density

rho = air_density(apogee)             # kg/m³
print(rho)

g = 9.81                  # m/s²
mass = 55.0               # rocket mass (kg)
m_p = 1.0                 # parachute mass (kg)
deploy_delay = 3.0        # s
t_reefed = 0.5           # partial inflation duration (s)
t_disreef = 1          # additional inflation duration (s)
dt = 0.001                # s

# Parachute and cord parameters
Cd_chute = 2.2
Cd_partial = 1.9
A_chute = ((3.0/2)**2)*(np.pi)
A_partial = ((2.0/2)**2)*(np.pi)  # effective area after reefed phase
k = 500.0                 # spring constant (N/m)
L0 = 10.0                 # unstretched cord length (m)

# Lists for results
time_list, rocket_alt_list, rocket_vel_list, rocket_acc_list = [], [], [], []
parachute_alt_list, parachute_vel_list, parachute_acc_list, tension_list = [], [], [], []

# Pre-deployment: rocket free-fall
t = 0.0
y_r = apogee
v_r = 0.0

while t < deploy_delay and y_r > 0:
    a_r = g
    time_list.append(t)
    rocket_alt_list.append(y_r)
    rocket_vel_list.append(v_r)
    rocket_acc_list.append(a_r)
    parachute_alt_list.append(np.nan)
    parachute_vel_list.append(np.nan)
    parachute_acc_list.append(np.nan)
    tension_list.append(0.0)
    
    v_r += a_r * dt
    y_r -= v_r * dt
    t += dt

# Deployment: initialize parachute state
y_p = y_r + L0
v_p = v_r

# Post-deployment: two-body simulation
while y_r > 2200:

    # Compute effective parachute area and drag coefficient with two-phase inflation
    if t < deploy_delay:
        A_eff = 0.0
    elif deploy_delay <= t < deploy_delay + t_reefed:
        A_eff = A_partial * ((t - deploy_delay) / t_reefed)
        Cd = Cd_partial
    elif deploy_delay + t_reefed <= t < deploy_delay + t_reefed + t_disreef:
        A_eff = A_partial + (A_chute - A_partial) * ((t - (deploy_delay + t_reefed)) / t_disreef)
        Cd = Cd_partial + (Cd_chute - Cd_partial) * ((t - (deploy_delay + t_reefed)) / t_disreef)
    else:
        A_eff = A_chute
        Cd = Cd_chute

    # Spring extension and tension
    d = y_p - y_r
    T = k * (d - L0) if d > L0 else 0.0

    # Forces and accelerations
    a_r = g - T / mass
    drag = 0.5 * rho * Cd * A_eff * (v_p**2)
    a_p = g + T / m_p - drag / m_p

    # Record data
    time_list.append(t)
    rocket_alt_list.append(y_r)
    rocket_vel_list.append(v_r)
    rocket_acc_list.append(a_r)
    parachute_alt_list.append(y_p)
    parachute_vel_list.append(v_p)
    parachute_acc_list.append(a_p)
    tension_list.append(T)

    # Update states
    v_r += a_r * dt
    y_r -= v_r * dt
    v_p += a_p * dt
    y_p -= v_p * dt
    t += dt

time_arr = np.array(time_list)
rocket_alt = np.array(rocket_alt_list)
rocket_vel = np.array(rocket_vel_list)
rocket_acc = np.array(rocket_acc_list)
parachute_alt = np.array(parachute_alt_list)
parachute_vel = np.array(parachute_vel_list)
parachute_acc = np.array(parachute_acc_list)
tension = np.array(tension_list)

plt.figure(figsize=(10,8))
plt.subplot(4, 1, 1)
plt.plot(time_arr, rocket_alt, label="Rocket")
plt.plot(time_arr, parachute_alt, label="Parachute", linestyle="--")
plt.ylabel("Altitude (m)")
plt.legend()

plt.subplot(4, 1, 2)
plt.plot(time_arr, rocket_vel, label="Rocket")
plt.plot(time_arr, parachute_vel, label="Parachute", linestyle="--")
plt.ylabel("Velocity (m/s)")
plt.legend()

plt.subplot(4, 1, 3)
plt.plot(time_arr, rocket_acc, label="Rocket")
plt.plot(time_arr, parachute_acc, label="Parachute", linestyle="--")
plt.ylabel("Acceleration (m/s²)")
plt.legend()

plt.subplot(4, 1, 4)
plt.plot(time_arr, tension)
plt.ylabel("Tension (N)")
plt.xlabel("Time (s)")
plt.tight_layout()
plt.show()

print("Maximum tension force:", np.max(tension), "N")
print("Maximum rocket deceleration:", np.max(rocket_acc), "m/s²")
