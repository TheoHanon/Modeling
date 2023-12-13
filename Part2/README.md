# Modelling & Analysis of Dynamical System : Project PART II

This project simulates tree cover across three different continents, aiming to model equilibrium, predict tree cover status by the year 2100, and examine the effects of various state variables. It also includes visualization tools for interactive analysis and phase portrait animations.

## Description

The questions are split into several Python scripts that perform different functions within the project:
- `q1.py`: Computes the equilibrium state of the tree cover model.
- `q2.py`: Uses historical data to make tree cover predictions for the year 2100.
- `q3.py`: Introduces additional state variables into the model and simulates their effects over time, including generating an animation.
- `main.py`: Serves as the entry point to launch the functions defined in `q1.py`, `q2.py`, and `q3.py`.
- `interactivePlotQ3.py`: Provides tools to interactively visualize the effect of parameters on the dynamics of the system (for question 3).
- `phasePortrait.py`: Generates phase portrait of our system (for question 3).

Data files are located within the `data` directory and are critical to the model's predictive functions.

## Getting Started

### Dependencies
Required libraries:
* matplotlib
* numpy
* scipy
* cartopy

### Executing Program

You can get the plots of the report and of the animation using:
```
python3 main.py
```
Caution the animation can take time  depending of your computer. If you want to interactively visualize the effect of parameters on the dynamics of the system covered in the video, you can launch it with
```
python3 interactivePlotQ3.py
```
and if you want to display the phase portrait of the new model (not covered in the video) for several initial conditions, you can launch it with 
```
python3 phasePortrait.py
```


