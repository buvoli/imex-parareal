# parareal
Matlab code for generating parareal stability, accuracy, and convergence regions and prototyping numerical experiments.

### directory structure

- `mains`: contains all the main scripts for generating the paper figures.
- `stability`: contains all the code for calculating stability and convergence properties of IMEX integrators and parareal.
- `stepper`: contains the code for running a *serial* parareal algorithm that is useful for testing combinations.

### Additional Stability and Convergence Plots
Additional stability plots are located in the directory `mains/grid-plots-stability/figures`.
Diagrams show all coarse/fine pairings using IMEX-RK methods of orders 1 to 4 at three different zoom levels.
