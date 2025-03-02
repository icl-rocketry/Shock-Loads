import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt

# --- Calculation Functions ---

def air_density(altitude):
    rho0 = 1.225  # Sea-level air density in kg/m^3
    T0 = 288.15  # Sea-level standard temperature in K
    Lapse_rate = 0.0065  # Temperature lapse rate in K/m
    R = 287.05  # Specific gas constant for dry air in J/(kg·K)
    g = 9.81  # m/s^2 (gravitational acceleration)
    
    T = T0 - Lapse_rate * altitude
    P = rho0 * np.exp(-g * altitude / (R * T))
    rho = P / (R * T)
    return max(rho, 0)

def terminal_velocity(mass, Cd, A, h=None, rho=None):
    g = 9.81  # m/s^2 (gravitational acceleration)
    if rho is None:
        if h is None:
            raise ValueError("Either 'h' or 'rho' must be provided")
        rho = air_density(h)
    return np.sqrt((2 * mass * g) / (rho * Cd * A))

def drag_force(v, Cd, A, rho=None, h=None):
    if rho is None:
        if h is None:
            raise ValueError("Either 'h' or 'rho' must be provided")
        rho = air_density(h)
    return 0.5 * rho * v**2 * Cd * A

def cord_stretch_damped(force, length, cross_sectional_area, modulus_of_elasticity, damping_coefficient, velocity):
    return (force * length) / (cross_sectional_area * modulus_of_elasticity) + damping_coefficient * velocity

# --- UI Calculation Function ---

def air_density_isa(altitude):
    # Constants for the International Standard Atmosphere (ISA)
    T0 = 288.15  # Sea-level standard temperature in K
    P0 = 101325  # Sea-level standard pressure in Pa
    L = 0.0065  # Temperature lapse rate in K/m
    R = 287.05  # Specific gas constant for dry air in J/(kg·K)
    g = 9.81  # m/s^2 (gravitational acceleration)
    
    if altitude < 11000:  # Troposphere
        T = T0 - L * altitude
        P = P0 * (T / T0) ** (g / (R * L))
    else:  # Above Troposphere (simplified)
        T = T0 - L * 11000
        P = P0 * (T / T0) ** (g / (R * L)) * np.exp(-g * (altitude - 11000) / (R * T))
    
    rho = P / (R * T)
    return max(rho, 0)

def calculate_and_display():
    try:
        # Retrieve input values from UI
        mass_initial = float(entry_mass.get())
        horizontal_speed = float(entry_horizontal_speed.get())
        parachute_diameter_full = float(entry_diameter_full.get())
        reefed_diameter = float(entry_reefed_diameter.get())
        drag_coefficient_full = float(entry_drag_full.get())
        reefed_drag_coefficient = float(entry_drag_reefed.get())
        reefed_duration = float(entry_reefed_duration.get())
        reefed_deployment_time = float(entry_deployment_time.get())  # time to deploy reefed chute
        full_deployment_time = float(entry_full_deployment_time.get())  # time to deploy full chute
        cord_length = float(entry_cord_length.get())
        cord_diameter = float(entry_cord_diameter.get())
        cord_yield_strength = float(entry_cord_yield_strength.get())
        cord_modulus_of_elasticity = float(entry_cord_modulus.get())
        safety_factor = float(entry_safety_factor.get())

        # Calculate areas from diameters
        parachute_area_full = np.pi * (parachute_diameter_full / 2) ** 2
        reefed_area = np.pi * (reefed_diameter / 2) ** 2

        # Constants
        g = 9.81
        cord_cross_sectional_area = np.pi * (cord_diameter / 2) ** 2

        # Simulation parameters
        time_step = 0.0001  # Time step for simulation in seconds
        total_time = 20  # Total simulation time in seconds

        # Initial conditions
        initial_velocity = 0  # m/s (initial velocity at apogee)
        initial_height = 9000  # m (initial height at apogee)
        mass = mass_initial  # Initial mass

        # Time array for simulation
        time_array = np.arange(0, total_time, time_step)

        # Arrays to store results
        velocity_array = np.zeros_like(time_array)
        height_array = np.zeros_like(time_array)
        drag_force_array = np.zeros_like(time_array)
        shock_load_array = np.zeros_like(time_array)
        mass_array = np.zeros_like(time_array)

        # Initial values
        velocity_array[0] = initial_velocity
        height_array[0] = initial_height
        mass_array[0] = mass

        # Simulation loop
        for i in range(1, len(time_array)):
            t = time_array[i]
            if t < reefed_duration:
                Cd = reefed_drag_coefficient
                A = reefed_area
            else:
                Cd = drag_coefficient_full
                A = parachute_area_full

            # Calculate height
            height_array[i] = height_array[i-1] - velocity_array[i-1] * time_step

            # Calculate air density at current height using ISA model
            rho = air_density_isa(height_array[i])

            # Calculate drag force
            drag_force_array[i] = drag_force(velocity_array[i-1], Cd, A, rho=rho)

            # Calculate acceleration
            acceleration = g - (drag_force_array[i] / mass)

            # Calculate velocity
            velocity_array[i] = velocity_array[i-1] + acceleration * time_step

            # Update mass (example: linear decrease over time)
            mass = mass_initial * (1 - t / total_time)
            mass_array[i] = mass

            # Calculate shock load on the shock cord
            shock_load_array[i] = cord_stretch_damped(drag_force_array[i], cord_length, cord_cross_sectional_area, cord_modulus_of_elasticity, 100, velocity_array[i])

        # Plot results
        plt.figure(figsize=(10, 6))
        plt.plot(time_array, shock_load_array, label='Shock Load on Shock Cord')
        plt.xlabel('Time (s)')
        plt.ylabel('Shock Load (N)')
        plt.title('Shock Load on Shock Cord During Parachute Deployment')
        plt.grid(True)
        plt.legend()
        plt.show()

    except ValueError:
        messagebox.showerror("Invalid input", "Please enter valid numeric values for all fields.")

# --- UI Setup ---

root = tk.Tk()
root.title("Shock Load Analysis")

# Input and Output Fields in Aligned Layout
inputs_outputs = [
    ("Rocket Mass (kg)", "55.0", None, None),
    ("Horizontal Speed at Apogee (m/s)", "0", None, None),
    ("Full Parachute Diameter (m)", "6.5", None, None),
    ("Reefed Parachute Diameter (m)", "1.6", None, None),
    ("Drag Coefficient (Full)", "2.2", None, None),
    ("Drag Coefficient (Reefed)", "2.0", None, None),
    ("Reefed Duration (s)", "5", None, None),
    ("Reefed Deployment Time (s)", "5", None, None),
    ("Full Deployment Time (s)", "2", None, None),
    ("Shock Cord Length (m)", "2.5", None, None),
    ("Shock Cord Diameter (m)", "0.01", None, None),
    ("Cord Yield Strength (Pa)", "5e8", None, None),
    ("Cord Elastic Modulus (Pa)", "4.7625e+03", None, None),
    ("Safety Factor", "2.0", None, None),
]

# Create input and output fields in aligned layout
entry_vars = {}
for i, (label_text, default_value, output_label, result_var) in enumerate(inputs_outputs):
    # Left Column: Input Label and Entry
    ttk.Label(root, text=label_text).grid(row=i, column=0, padx=5, pady=2, sticky="w")
    entry = ttk.Entry(root)
    entry.insert(0, default_value)
    entry.grid(row=i, column=1, padx=5, pady=2)
    entry_vars[label_text] = entry

# Reference input fields for calculation
entry_mass = entry_vars["Rocket Mass (kg)"]
entry_horizontal_speed = entry_vars["Horizontal Speed at Apogee (m/s)"]
entry_diameter_full = entry_vars["Full Parachute Diameter (m)"]
entry_reefed_diameter = entry_vars["Reefed Parachute Diameter (m)"]
entry_drag_full = entry_vars["Drag Coefficient (Full)"]
entry_drag_reefed = entry_vars["Drag Coefficient (Reefed)"]
entry_reefed_duration = entry_vars["Reefed Duration (s)"]
entry_deployment_time = entry_vars["Reefed Deployment Time (s)"]
entry_full_deployment_time = entry_vars["Full Deployment Time (s)"]
entry_cord_length = entry_vars["Shock Cord Length (m)"]
entry_cord_diameter = entry_vars["Shock Cord Diameter (m)"]
entry_cord_yield_strength = entry_vars["Cord Yield Strength (Pa)"]
entry_cord_modulus = entry_vars["Cord Elastic Modulus (Pa)"]
entry_safety_factor = entry_vars["Safety Factor"]

# Calculate Button
calculate_button = ttk.Button(root, text="Calculate", command=calculate_and_display)
calculate_button.grid(row=len(inputs_outputs), column=1, pady=10)

root.mainloop()

