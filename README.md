# ultrasound-simulations
Simplified interfaces for ultrasound simulations

All the codes in this repo are ment to provide an OOP interface to create, save and visualize ultrasound simulations.

## The available models
### Field2Simulation
An interface for the [field-II](http://field-ii.dk//) matlab package, which finds the pressure field profile in a given set of coordinates
by solving the Fersnel integral (wave propagation model).

### kWave_linear_array_in_medium
An interface for the [k-Wave](http://www.k-wave.org/) matlab toolbox which solves the wave propagation with tracking of particle speed.
* This model provides a more smooth solution for the pressure field but requires more computation time compared to Field-II

### Marmottant_model
An interface for Marmottant model for bubble oscillation, buckling and rupture in response to oscillatory prssure changes.
[read more about Marmottant model](https://asa.scitation.org/doi/10.1121/1.2109427)

* Important! make sure to download the required packages and add them to your MATLAB path.
