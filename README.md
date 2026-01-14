# Seasonal-Optimal-Control
This repository provides the MATLAB implementation of the climate-driven malaria model with seasonal optimal control exactly as used in the manuscript. The code is split into modular files for clarity: baseline simulation, optimal control (forward–backward sweep), climate forcing, and plotting. Running main.m reproduces all figures and results in the paper.


The provided code implements three time-dependent controls:

u1 : personal protection

u2 : treatment

u3 : vector control

By default, the script runs the full optimal control strategy where all three controls are active.





No other changes to the code are needed.


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
