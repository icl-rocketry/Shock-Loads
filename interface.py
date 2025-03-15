import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt

def air_density(h):
    t0 = 288.15
    pho0 = 1.225
    r = 287.058
    P = 101325
    lapse_rate = 0.0065
    g = 9.81

    T = t0 - lapse_rate * h
    pressure = P * (T / t0) ** (g / (lapse_rate * r))
    density = pressure / (r * T)
    return density

def simulate():
    try:
        apogee = float(entry_apogee.get())
        mass = float(entry_mass.get())
        m_p = float(entry_parachute_mass.get())
        deploy_delay = float(entry_deploy_delay.get())
        t_reefed = float(entry_t_reefed.get())
        t_disreef = float(entry_t_disreef.get())
        dt = float(entry_dt.get())
        Cd_chute = float(entry_Cd_chute.get())
        Cd_partial = float(entry_Cd_partial.get())
        A_chute = np.pi * (float(entry_diameter_chute.get()) / 2) ** 2
        A_partial = np.pi * (float(entry_diameter_partial.get()) / 2) ** 2
        k = float(entry_k.get())
        L0 = float(entry_L0.get())

        rho = air_density(apogee)
        g = 9.81

        time_list, rocket_alt_list, rocket_vel_list, rocket_acc_list = [], [], [], []
        parachute_alt_list, parachute_vel_list, parachute_acc_list, tension_list = [], [], [], []

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

        y_p = y_r + L0
        v_p = v_r

        while y_r > 2200:
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

            d = y_p - y_r
            T = k * (d - L0) if d > L0 else 0.0

            a_r = g - T / mass
            drag = 0.5 * rho * Cd * A_eff * (v_p**2)
            a_p = g + T / m_p - drag / m_p

            time_list.append(t)
            rocket_alt_list.append(y_r)
            rocket_vel_list.append(v_r)
            rocket_acc_list.append(a_r)
            parachute_alt_list.append(y_p)
            parachute_vel_list.append(v_p)
            parachute_acc_list.append(a_p)
            tension_list.append(T)

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

        plt.figure(figsize=(10, 8))

        ax1 = plt.subplot(4, 1, 1)
        ax1.plot(time_arr, rocket_alt, color='xkcd:dark blue', label="Rocket")
        ax1.plot(time_arr, parachute_alt, color='xkcd:red', label="Parachute", linestyle=":")
        ax1.set_ylabel("Altitude (m)")
        ax1.legend()
        ax1.grid(which='major', color='xkcd:dark blue', alpha=0.2, linewidth=0.4)
        ax1.grid(which='minor', color='xkcd:dark blue', alpha=0.1, linewidth=0.2)
        ax1.minorticks_on()
        ax1.tick_params(which='both', direction='in', top=True, right=True)

        ax2 = plt.subplot(4, 1, 2)
        ax2.plot(time_arr, rocket_vel, color='xkcd:dark blue', label="Rocket")
        ax2.plot(time_arr, parachute_vel, color='xkcd:red', label="Parachute", linestyle=":")
        ax2.set_ylabel("Velocity (m/s)")
        ax2.legend()
        ax2.grid(which='major', color='xkcd:dark blue', alpha=0.2, linewidth=0.4)
        ax2.grid(which='minor', color='xkcd:dark blue', alpha=0.1, linewidth=0.2)
        ax2.minorticks_on()
        ax2.tick_params(which='both', direction='in', top=True, right=True)

        ax3 = plt.subplot(4, 1, 3)
        ax3.plot(time_arr, rocket_acc, color='xkcd:dark blue', label="Rocket")
        ax3.plot(time_arr, parachute_acc, color='xkcd:red', label="Parachute", linestyle=":")
        ax3.set_ylabel("Acceleration (m/s²)")
        ax3.legend()
        ax3.grid(which='major', color='xkcd:dark blue', alpha=0.2, linewidth=0.4)
        ax3.grid(which='minor', color='xkcd:dark blue', alpha=0.1, linewidth=0.2)
        ax3.minorticks_on()
        ax3.tick_params(which='both', direction='in', top=True, right=True)

        ax4 = plt.subplot(4, 1, 4)
        ax4.plot(time_arr, tension, color='xkcd:dark blue')
        ax4.set_ylabel("Tension (N)")
        ax4.set_xlabel("Time (s)")
        ax4.grid(which='major', color='xkcd:dark blue', alpha=0.2, linewidth=0.4)
        ax4.grid(which='minor', color='xkcd:dark blue', alpha=0.1, linewidth=0.2)
        ax4.minorticks_on()
        ax4.tick_params(which='both', direction='in', top=True, right=True)

        plt.tight_layout()
        plt.show()

        messagebox.showinfo("Simulation Complete", f"Maximum tension force: {np.max(tension):.2f} N") #\nFor safety factor of {safety_factor}: {np.max(tension_sf):.2f} N")

    except ValueError:
        messagebox.showerror("Invalid input", "Please enter valid numeric values for all fields.")

def set_nimbus24_values():

    entry_apogee.delete(0, tk.END)
    entry_apogee.insert(0, "3000")

    entry_mass.delete(0, tk.END)
    entry_mass.insert(0, "55")

    entry_parachute_mass.delete(0, tk.END)
    entry_parachute_mass.insert(0, "1.0")

    entry_deploy_delay.delete(0, tk.END)
    entry_deploy_delay.insert(0, "10")

    entry_t_reefed.delete(0, tk.END)
    entry_t_reefed.insert(0, "0.5")

    entry_t_disreef.delete(0, tk.END)
    entry_t_disreef.insert(0, "1.0")

    entry_Cd_chute.delete(0, tk.END)
    entry_Cd_chute.insert(0, "2.2") # didn't disreef

    entry_Cd_partial.delete(0, tk.END)
    entry_Cd_partial.insert(0, "2.2")

    entry_diameter_chute.delete(0, tk.END)
    entry_diameter_chute.insert(0, "2.0") # didn't disreef

    entry_diameter_partial.delete(0, tk.END)
    entry_diameter_partial.insert(0, "2.0")

    entry_k.delete(0, tk.END)
    entry_k.insert(0, "500")

    entry_L0.delete(0, tk.END)
    entry_L0.insert(0, "10") # deployed late


root = tk.Tk()
root.title("Parachute Deployment Simulation")

inputs = [
    ("Apogee (m)", "3000.0"),
    ("Rocket Mass (kg)", "55.0"),
    ("Parachute Mass (kg)", "1.0"),
    ("Deployment Delay (s)", "3.0"),
    ("Reefed Inflation Duration (s)", "0.5"),
    ("Disreef Inflation Duration (s)", "1.0"),
    ("Time Step (s)", "0.001"),
    ("Drag Coefficient (Full)", "2.2"),
    ("Drag Coefficient (Partial)", "1.9"),
    ("Parachute Diameter (m)", "3.0"),
    ("Partial Parachute Diameter (m)", "2.0"),
    ("Spring Constant (N/m)", "500.0"),
    ("Unstretched Cord Length (m)", "10.0"),
]

entry_vars = {}
for i, (label_text, default_value) in enumerate(inputs):
    ttk.Label(root, text=label_text).grid(row=i, column=0, padx=5, pady=2, sticky="w")
    entry = ttk.Entry(root)
    entry.insert(0, default_value)
    entry.grid(row=i, column=1, padx=5, pady=2)
    entry_vars[label_text] = entry

entry_apogee = entry_vars["Apogee (m)"]
entry_mass = entry_vars["Rocket Mass (kg)"]
entry_parachute_mass = entry_vars["Parachute Mass (kg)"]
entry_deploy_delay = entry_vars["Deployment Delay (s)"]
entry_t_reefed = entry_vars["Reefed Inflation Duration (s)"]
entry_t_disreef = entry_vars["Disreef Inflation Duration (s)"]
entry_dt = entry_vars["Time Step (s)"]
entry_Cd_chute = entry_vars["Drag Coefficient (Full)"]
entry_Cd_partial = entry_vars["Drag Coefficient (Partial)"]
entry_diameter_chute = entry_vars["Parachute Diameter (m)"]
entry_diameter_partial = entry_vars["Partial Parachute Diameter (m)"]
entry_k = entry_vars["Spring Constant (N/m)"]
entry_L0 = entry_vars["Unstretched Cord Length (m)"]

calculate_button = ttk.Button(root, text="Simulate", command=simulate)
calculate_button.grid(row=len(inputs), column=1, pady=10)

nimbus24_button = ttk.Button(root, text="Nimbus24", command=set_nimbus24_values)
nimbus24_button.grid(row=len(inputs)+1, column=1, pady=10)

root.mainloop()