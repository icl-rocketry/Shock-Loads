import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
time_step = 0.0001  # Time step for simulation in seconds
total_time = 20  # Total simulation time in seconds

# Constants
m = 55  # kg (rocket mass)
g = 9.81  # m/s^2 (gravitational acceleration)
h = 9000  # m (apogee height)

# Drag coefficients
Cd_reefed = 2.0
Cd_deployed = 2.2

# Parachute diameters and areas
d_reefed = 1.6  # m
d_deployed = 6.5  # m
A_reefed = np.pi * (d_reefed / 2) ** 2
A_deployed = np.pi * (d_deployed / 2) ** 2

t_disreef = 5  # Time at which the parachute is fully deployed in seconds (t=0 is apogee)

# Shock cord parameters
L = 2.5  # m (tether length)
k = 4.7625e+03  # N/m (spring stiffness)
c = 100  # Ns/m (damping coefficient)


## Density calcs

# Air properties
rho0 = 1.225  # Sea-level air density in kg/m^3
T0 = 288.15  # Sea-level standard temperature in K
Lapse_rate = 0.0065  # Temperature lapse rate in K/m
R = 287.05  # Specific gas constant for dry air in J/(kg·K)

# ISA model for air density as a function of altitude
def air_density(altitude):
    T = T0 - Lapse_rate * altitude
    P = rho0 * np.exp(-g * altitude / (R * T))
    rho = P / (R * T)
    return max(rho, 0)


# Calcs

def terminal_velocity(mass, Cd, A, h=None, rho=None):
    if rho is None:
        if h is None:
            raise ValueError("Either 'h' or 'rho' must be provided")
        rho = air_density(h)
    return np.sqrt((2 * mass * g) / (rho * Cd * A))

def deployment_velocity(t, h_speed):
    """t=0 is apogee"""
    return np.sqrt(h_speed**2 + (g * t)**2)

def drag_force(v, Cd, A, rho=None, h=None):
    if rho is None:
        if h is None:
            raise ValueError("Either 'h' or 'rho' must be provided")
        rho = air_density(h)
    return 0.5 * rho * v**2 * Cd * A

def cord_stretch(force, length, cross_sectional_area, modulus_of_elasticity):
    return (force * length) / (cross_sectional_area * modulus_of_elasticity)

def cord_stretch_damped(force, length, cross_sectional_area, modulus_of_elasticity, damping_coefficient, velocity):
    return (force * length) / (cross_sectional_area * modulus_of_elasticity) + damping_coefficient * velocity

def cord_stress(force, cross_sectional_area):
    return force / cross_sectional_area

def velocity_at_time(mass, Cd, A, t, h_speed, h=None, rho=None):
    if rho is None:
        if h is None:
            raise ValueError("Either 'h' or 'rho' must be provided")
        rho = air_density(h)
    term_velocity = terminal_velocity(mass, Cd, A, rho)
    tau = (2 * mass) / (rho * Cd * A * term_velocity)
    return h_speed * np.exp(-t / tau) + term_velocity * (1 - np.exp(-t / tau))


# Simulation

# Initial conditions
initial_velocity = 0  # m/s (initial velocity at apogee)
initial_height = h  # m (initial height at apogee)

# Time array for simulation
time_array = np.arange(0, total_time, time_step)

# Arrays to store results
velocity_array = np.zeros_like(time_array)
height_array = np.zeros_like(time_array)
drag_force_array = np.zeros_like(time_array)
shock_load_array = np.zeros_like(time_array)

# Initial values
velocity_array[0] = initial_velocity
height_array[0] = initial_height

# Simulation loop
for i in range(1, len(time_array)):
    t = time_array[i]
    if t < t_disreef:
        Cd = Cd_reefed
        A = A_reefed
    else:
        Cd = Cd_deployed
        A = A_deployed

    # Calculate height
    height_array[i] = height_array[i-1] - velocity_array[i-1] * time_step

    # Calculate air density at current height
    rho = air_density(height_array[i])

    # Calculate drag force
    drag_force_array[i] = drag_force(velocity_array[i-1], Cd, A, rho=rho)

    # Calculate acceleration
    acceleration = g - (drag_force_array[i] / m)

    # Calculate velocity
    velocity_array[i] = velocity_array[i-1] + acceleration * time_step

    # Calculate shock load on the shock cord
    shock_load_array[i] = cord_stretch_damped(drag_force_array[i], L, np.pi * (d_reefed / 2)**2, k, c, velocity_array[i])

# Plot results
plt.plot(time_array, shock_load_array, label='Shock Load on Shock Cord')
plt.xlabel('Time (s)')
plt.ylabel('Shock Load (N)')
plt.title('Shock Load on Shock Cord During Parachute Deployment')
plt.grid(True)
plt.show()

