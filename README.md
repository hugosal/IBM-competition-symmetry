# IBM-competition-symmetry
Codes and data analysis related to "What can plant size variation tell us about symmetry of competition?"
Salinas Hugo, Veneklaas Erik J., Renton Michael 2023.


## Overview
This is a model that simulates the development of a population of plants to study. 
It is used to study how to identify competition symmetry from plant size variation.


### Finding configurations of points
One of the requirements of the model is to have coordinates of N points in 2D torus square space such that every point is at the same distance from the closest 4 points.
To find configurations of points that meet this requirement we used an iterative process of applying roations and scaling to a set of points in a grid.
To generate the *point_configurations.csv* file with the coordinates of the points necessary to run simulations of plant populations:

```bash
python3 looping_finding_points_config.py
mkdir configurations
mv point_configurations.csv configurations/
```

### Runs of simulations
The R package *CirclesIntersections*, and the *auxiliary_functions.R* script are needed.
The conditions and number in which the simulations will be made are hard coded in the Rscript but can be modified.


### Data analysis

*spatial_pattern_randomness_tests.R* Analyse competition experiments, PENDING TO DOCUMENT
*make_development_animation.R* efefefefefef
*ibm facilitation_output_analysis.R* wefefefefeff


## Authors

Elizabeth Trevenen and Michael Renton developed the model. Contributions to the model, experiments and data analysis made by Hugo Salinas.


## License

[MIT](https://github.com/hugosal/IBM-competition-symmetry/blob/main/LICENSE)

