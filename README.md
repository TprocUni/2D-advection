

### In-Depth Overview (

This program simulates a **2D advection** process, in which a Gaussian distribution of values is moved or "advected" across a grid in the x and y directions based on specified velocity fields. Advection, in a broader sense, refers to the transport of a substance or property through a fluid or field due to the motion of the fluid itself. This simulation can be useful in understanding how particles or substances disperse in various fields, such as fluid dynamics or atmospheric sciences.

#### Initialization and Configuration

1. **Grid Properties**:
   - The grid size is defined by the constants `NX` and `NY`, representing the number of points in the x and y directions, respectively.
   - The grid spans specific x and y boundaries, with `xmin`, `xmax`, `ymin`, and `ymax` defining the minimum and maximum extents of the simulation space.
   - The program calculates `dx` and `dy`, which are the distances between consecutive grid points in x and y, respectively. These values are crucial for calculating changes in the field over time, allowing the program to satisfy the **Courant–Friedrichs–Lewy (CFL) condition** for stability.

2. **Gaussian Initial Conditions**:
   - The initial state of the grid, `u`, is set up as a **Gaussian distribution** centered at `(x0, y0)` with specified widths `sigmax` and `sigmay` for x and y, respectively.
   - A Gaussian profile is chosen for its smooth, bell-shaped properties, which make it ideal for observing how advection alters the distribution over time.

3. **Boundary Conditions**:
   - Boundary conditions are set on all four edges of the grid to maintain stability during simulation.
   - The boundaries (`bval_left`, `bval_right`, `bval_lower`, and `bval_upper`) are given a constant value of zero, ensuring that values at the grid edges remain fixed and do not influence the interior values as they advect.

#### Velocity Field and CFL Condition

- The **velocity field** for advection is set using `velx` and `vely`. The x-direction velocity, `velx`, varies with y, following a **shear velocity profile** calculated by the function `calc_x_shear()`.
- The function `calc_x_shear()` uses a **logarithmic profile** based on friction and roughness length to simulate horizontal shear, where velocity increases with height.
- **CFL Condition**: The time step `dt` is determined by the CFL condition, which ensures numerical stability by limiting the time step based on grid spacing and velocity. This condition is critical to prevent the simulation from becoming unstable due to high velocity or small grid spacing.

#### Time Stepping and Parallelization

1. **Main Advection Loop**:
   - The program iterates over a specified number of time steps (`nsteps`), updating the values of `u` in each iteration to simulate the movement of the Gaussian distribution.
   - For each time step:
     - **Boundary Conditions** are applied to ensure stability.
     - The **rate of change** of `u` (stored in `dudt`) is calculated using **finite differences**. This approach involves subtracting neighboring values to determine the change at each point, approximating the partial derivative.
     - Finally, the values of `u` are updated based on `dudt` and `dt`.

2. **Parallelization**:
   - The program uses **OpenMP** to parallelize several loops. OpenMP is a library for parallel processing in C, which allows code to run across multiple CPU cores.
   - Loops that initialize the grid points (`x` and `y`), set up the initial Gaussian, and calculate `dudt` are parallelized, improving performance on systems with multiple cores.
   - However, some loops, such as those writing to files, cannot be parallelized due to **output dependencies**. Writing in parallel could cause data to be written out of order, resulting in incorrect or jumbled outputs.

#### Output Files

The program produces three output files:

1. **`initial.dat`**: Stores the initial values of `u` across the grid. This file provides a reference for how the Gaussian profile looks at the start of the simulation.
2. **`final.dat`**: Contains the values of `u` after all time steps, showing the final state of the advection process.
3. **`vertically_averaged.dat`**: This file provides a vertically averaged distribution of `u` across each x-coordinate. It averages `u` values along each column to show the mean distribution along the x-axis, giving a summary of how the Gaussian profile has dispersed vertically.

#### Performance Monitoring

The program records **elapsed time** using OpenMP's `omp_get_wtime()`, which marks the start and end times of the main computation. This time measurement helps monitor performance, especially useful for evaluating how parallelization affects the runtime.

#### Summary

In summary, this C program simulates the 2D advection of a Gaussian distribution, using finite differences and parallelization for efficiency. The shear velocity profile,
