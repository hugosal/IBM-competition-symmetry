# IBM-competition-symmetry
Codes and data analysis related to "What can plant size variation tell us about symmetry of competition?"
Salinas Hugo, Veneklaas Erik J., Renton Michael 2023.


## Overview
This is a model that simulates the development of a population of plants to study. 
It is used to study how to identify competition symmetry from plant size variation.


### Finding configurations of points
One of the requirements of the model is to have coordinates of N points in 2D torus square space such that every point is at the same distance from the closest 4 points.
To find configurations of points that meet this requirement we used an iterative process of applying roations and scaling to a set of points in a grid.
The generate the .csv file with the coordinates of the points necessary to run simulations of plant populations:

```bash
python3 looping_finding_points_config.py
```


### Data analysis

*analyse_community_experiments.R* Analyse competition experiments, PENDING TO DOCUMENT

*analyze_genotypes_over_time.R* Can be used to plot and test on biomass, and morphology metrics over time, using plants growing alone. 

*competitoin_experimentsanalyzer_homevs_visitor.R* Heatmaps and line plots of shoot biomass depending on the number of neighbours.

*genotype_library_summary.R* Basic inspection of root architecture parameters in a genotype pool.


*parameter_analyzer_trhough_time.R* Multidimensional analysis of the root architecture parameters as they change over time



## Authors

Elizabeth Trevenen and Michael Renton developed the model. Contributions to the model, experiments and data analysis made by Hugo Salinas.


## License

[MIT](https://github.com/hugosal/IBM-competition-symmetry/blob/main/LICENSE)

