# JMDT: Jank Mission Design Tool.
Exercise caution and a healthy dose of scientific skepticism.

# Features
- C++ code with a Python wrapper for plots and glue. It's very fast: about 800 simulated days in orbit per second using the most simple settings (and a time step of 10 seconds), and a respectable 150 simulated days per second using WGS84, and simulating both solar panels and drag.
- Takes 3D models of an arbitrary satellite as input, so it can figure out what power the satellite is getting from any orientation, and what the frontal area for drag is. Painfully enough, it even takes into account when parts of the satellite are blocking part of the solar panels. This is computed on the GPU using OpenCL.
- Three gravity models: point mass, WGS84, EGM96.
- Three atmospheric models: none, US1976 and NRLMSISE-00.
- Adams-Bashforth-Moulton integrator.
- Can simulate two satellites at a time.
- Change the control "mode": pointing at the other satellite, point radially, point prograde, point at the Sun. This changes the solar output and total drag, which can be useful when designing a mission and evaluating the power requirements and availability.
- When the time comes, one can *somewhat* trivially perform hardware in the loop with this tool as the physical plant and, before that, debug control algorithms.

# Contact
Joan Creus-Costa `<jcreus@stanford.edu>`
Sasha Maldonado `<amaldona@stanford.edu>`, but he's hardly as useful
