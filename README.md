# explicitTraps

-- Synopsis --
This is a program developed to simulate the transport, decay of, and pairwise interaction between excitons in two-dimensional materials. The simulation was developed for excitons in monolayer MoS_2, a 2D semiconductor. Excitons move diffusively but are beleagured by the presence of stationary deep traps. They trap and thermally detrap. When two excitons encounter each other, one is annihilated dissipating its energy as heat. Code in this repository outputs the time-dependent exciton population (as a function of density) in a transient spectroscopy experiment, and also the power-dependent steady-state quantum yield. 


-- Code Example --
newDiffSim.cpp - 
This program runs a simulation of excitons diffusing in a 2D semiconductor. The excitons are excited by a pulsed laser with a repetition rate and spot size chosen by the user. The simulation outputs the time (measured with respect to the last laser pulse) and spatial position of each exciton radiative decay event. This allows for the quantitative simulation of a time- and spatially-resolved photoluminescence experiment (see for example doi.org/10.1021/nl501190s). The simulation also captures exciton-exciton annihliation interactions and trapping-detrapping events.

QY.cpp - 
This program runs a simulation of excitons diffusing in a 2D semiconductor. The excitons are excited by a CW laser at random with an average excitation rate R. The simulation tracks when excitons decay radiatively and how many were excited allowing the calculation of the radiative quantum yield. The program simulates multiple excitation rates R to show the deleterious effects of exciton-exciton annihilation on the exciton quantum yield. This allows one to simulate steady-state PL quantum yield experiments (see for example doi.org/10.1126/science.aad2114)


-- Motivation --
2D semiconductors such as transition metal dichalcogenides are promising candidates for next-generation photodetectors, light emitting devices, and transistors. Excitons move diffusively in these materials. The energy landscape seen by a diffusing exciton can be riddled with deep traps, which complicate the diffusive transport such that traditional analytical diffusion equations don't capture the nuances of exciton transport (see doi.org/10.1103/PhysRevB.96.121404). Furthermore, exciton-exciton interactions limit quantum yields and consequently device operating efficiencies at modest exciton densities. Understanding exciton transport and interaction is important to designing better devices.

These simulations capture the nuances and complications of exciton transport and annihilation in these technologically compelling materials systems.


-- contact __
These are simulations for a relatively niche research field that I'm publishing to show my work. If you're interested in fixing/adpating this code, please reach out at aarongoodm@gmail.com

 
