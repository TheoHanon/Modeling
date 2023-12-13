# Modelling & Analysis of Dynamical System : Project PART II

This project simulates tree cover across three different continents, aiming to model equilibrium, predict tree cover status by the year 2100, and examine the effects of various state variables. It also includes visualization tools for interactive analysis and phase portrait animations.

## Description

The simulation is split into several Python scripts that perform different functions within the project:
- `q1.py`: Computes the equilibrium state of the tree cover model.
- `q2.py`: Uses historical data to make tree cover predictions for the year 2100.
- `q3.py`: Introduces additional state variables into the model and simulates their effects over time, including generating an animation.
- `main.py`: Serves as the entry point to launch the functions defined in `q1.py`, `q2.py`, and `q3.py`.
- `interactivePlotQ3.py`: Provides tools for interactive analysis of the simulation.
- `phasePortrait.py`: Generates phase portrait animations based on the simulation data.

Data files are located within the `data` directory and are critical to the model's predictive functions. All the plot of the report are located in `plot` directory. The report and the video are in the `report & video` directory.

## Getting Started

### Dependencies
Required libraries:
*  matplotlib
* numpy
* scipy
* cartopy

### Executing Program

You can execute the program in your terminal using:
```
python3 main.py
```

If you want to run the iteractive plot of our new model (covered in the video) you can launch it with
```
python3 interactivePlotQ3.py
```
and if you want to display some phase portrait of the new model (not covered in the video) you can launch it with 
```
python3 phasePortrait.py
```


