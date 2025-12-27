# Seasonal-Optimal-Control
This repository provides the MATLAB implementation of the climate-driven malaria model with seasonal optimal control exactly as used in the manuscript. The code is split into modular files for clarity: baseline simulation, optimal control (forward–backward sweep), climate forcing, and plotting. Running main.m reproduces all figures and results in the paper.


The provided code implements three time-dependent controls:

u1 : personal protection

u2 : treatment

u3 : vector control

By default, the script runs the full optimal control strategy where all three controls are active.

You can also simulate additional scenarios such as:

Only u1 active (u2 = 0, u3 = 0)

Only u2 active (u1 = 0, u3 = 0)

Only u3 active (u1 = 0, u2 = 0)

Any two controls active (for example, u1 and u2 active, u3 = 0)

How to Modify the Code

Inside the forward–backward sweep loop, after the line where the updated controls u1(i), u2(i), and u3(i) are computed, you can force the inactive controls to zero.

Examples:

To run the u1-only scenario:

u2(i) = 0;
u3(i) = 0;


To run the u2-only scenario:

u1(i) = 0;
u3(i) = 0;


To run the u3-only scenario:

u1(i) = 0;
u2(i) = 0;


To run a two-control scenario (example: u1 and u2 active):

u3(i) = 0;


No other changes to the code are needed.

Plotting

After modifying the controls as shown above, simply run the script again.
The code will automatically generate updated plots for:

susceptible, exposed, infected, and recovered hosts

aquatic, susceptible, and infected vectors

the time profiles of u1, u2, and u3

This allows users to reproduce all alternative control strategies used in the study without requiring additional scripts.


Model Summary

The model includes human (Sh, Eh, Ih, Rh) and mosquito (Am, Sm, Im) compartments with temperature-, rainfall-, and humidity-dependent parameters.
Optimal control uses the forward–backward sweep method with three controls: personal protection, treatment, and vector control.

Citation:If you use this code, please cite the associated manuscript.
Jain, H., Raidas, S. K., & Sinha, A. K. (2026). A climate-based malaria transmission model with seasonal optimal control and cost-effective analysis. Chaos, Solitons & Fractals, 204, 117802. https://doi.org/10.1016/j.chaos.2025.117802

Contact

Himanshu Jain
himanshujain.nitrr@gmail.com
GitHub: https://github.com/itsjain
