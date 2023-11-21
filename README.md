# BISCUIT

<img src="Figures/Biscuit-logo.jpeg" width="150" align="right">

BISCUIT (Bayesian Integration of Sequences Using Chronology Intersection Tie-points) is an R package for integrating the chronology of an input paleoclimate record with a target record that has known chronological control points.

It allows users to interactively select matched tie-points between the input and target records. BISCUIT then runs a Bayesian model to estimate the age-depth relationship of the input record by incorporating information from:

- The posterior distributions of ages at the tie-point depths in the target record
- Prior information about sediment accumulation rates and memory

This results in a full probabilistic age model for the input record calibrated to the chronology of the target.

***

## Installation

### Requirements

BISCUIT requires:

- R (version 3.6 or higher) 
- The following files from this GitHub repo:
  - biscuit.R
  - twalk.R

To install, clone or download the repository so the files are in your working directory.

***

## Usage 

1. **Prepare input data**

   Create a folder `input-target` containing:
   
   - `input_proxy.csv` - Input proxy data 
   - `target_proxy.csv` - Target proxy data
   - `target/` - Output folder from [Bacon](https://github.com/Maarten14C/rbacon.git) or [Plum](https://github.com/Maarten14C/rplum.git)

2. **Run BISCUIT**

   ```R
   output <- Biscuit(Input="input", Target="target", 
                    folder='data/',  
                    n_tie_points=5,
                    run_target=FALSE)
   ```
   
3. **Interactively select tie points** between the input and target records

4. **Output** age model, tie points, and diagnostics are returned for the input record  
   
See [examples](examples) for sample data and usage.

*** 

## Function Arguments

`Biscuit()` accepts the following arguments:

- `Input` (str): Name of input record CSV file
- `Target` (str): Name of target record CSV file 
- `folder` (str): Path to `input-target` folder
- `n_sections`: Number of depth sections to model. If TRUE, sections are automatically calculated.

- `n_tie_points`: Number of tie points to select between input and target records.

- `sampling_date`: Date the input record was collected. 

- `thin`: Thinning interval for MCMC sampling.

- `burn`: Number of burn-in iterations for MCMC sampling. 

- `iters`: Total number of MCMC iterations.

- `shape_acc`: Shape parameter for prior on accumulation rate.

- `mean_acc`: Mean parameter for prior on accumulation rate.

- `strength_mem`: Strength hyperparameter for prior on sediment memory.

- `mean_mem`: Mean hyperparameter for prior on sediment memory.

- `run_target`: Whether to run Bacon/Plum on target chronology.



***

## Output 

`Biscuit()` returns a list containing:

- `age_models`: Matrix of age models at depth sections
- `elbows`: Depth sections
- `age_means`: Mean age at depth sections 
- `tie_points`: Dataframe of selected tie points
- `age_min`: 5th percentile age at depths
- `age_medium`: 50th percentile (median) age at depths  
- `age_max`: 95th percentile age at depths

***

<!-- ## Examples

See the [examples/](examples) folder for sample data and scripts. -->

## Authors

BISCUIT was developed by [Marco A. Aquino-Lopez](aquino@cimat.mx).

## Issues

Please submit bugs and feature requests via GitHub issues or by emailing the [author](aquino@cimat.mx).

## License 

BISCUIT is licensed under the MIT License. 

## Acknowledgements

We acknowledge contributors who assisted with development of BISCUIT.