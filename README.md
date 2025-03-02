# Shock Load Analysis

Application to simulate the shock load on a parachute's shock cord during deployment, considering varying air densities at different altitudes and dynamic system properties.

## Requirements

- Python
- NumPy
- Matplotlib
- Tkinter (for user interface)


## Usage

1. Run the script:
   ```bash
   python interface.py
   ```

2. Input Parameters

3. Click the **Calculate** button to start the simulation.

### Notes
- The simulation uses the International Standard Atmosphere (ISA) model to calculate air density.
- The shock cord stretch is modelled using Hooke's Law with additional damping.