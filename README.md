# bbh-subpopulations-scripts

Scripts for running the binned gaussian process analyses presented in [2404.03166](https://arxiv.org/abs/2404.03166). To run these scripts the following packages are necessary:

- https://github.com/AnaryaRay1/bbh-subpopulations-scripts/tree/main ([9f853ce](https://github.com/AnaryaRay1/gppop/commit/9f853cecc3ab1c5ee4358baaf6876492d59c1945))
- https://git.ligo.org/reed.essick/gw-distributions.git ([5d3804](https://git.ligo.org/reed.essick/gw-distributions/-/commit/5d38045f93cf5367ff031689ea27e1d640f0f3c1))
- https://pypi.org/project/bilby/ (latest)
- https://pypi.org/project/gwpopulation/ (latest)

as well as the datasets uploaded to: https://zenodo.org/records/12746314. Please see the following instructions (for a linux-based os) for reproducing the results of the paper.

### Set up environment

After creating and activating a python virtual environment, and clonning this repo take the following steps:

```
$ git clone -b spin-dev https://github.com/AnaryaRay1/gppop.git
$ cd gppop
$ pip install .
$ cd ..
$ git clone https://git.ligo.org/reed.essick/gw-distributions.git
$ cd gw-distributions
$ git apply /path/to/bbh-subpopulations-scripts/gwdistributions.patch --quiet
$ pip install .
$ cd ..
$ pip install bilby
$ pip install gwpopulations
```

**Note**, the reason we need to patch ```gw-distributions``` is as follows. We use this package to calculate the marginal prior on effective spins corresponding to the full spin priors used in parameter estimation and in drawing the sensitivity injections. For GWTC-3, this prior corresponds to distribution of component spins that is isotropic, and unifrom in margnitude, which is derived in Eq 10 of https://arxiv.org/abs/2104.09508 and coded up in the master branch of ```gw-distributions```. However, for the simulated events, we use an aligned prior on component spins, i.e. Eq. 7 of https://arxiv.org/abs/2104.09508, which is not coded up in master and hence needs to be patched into the local copy.

### Download data 

To download the required data from the data-release:

```
$ mkdir -p /path/to/data-dir/
$ cd /path/to/bbh-subpopulations-scripts/
$ ./download_data /path/to/data-dir/
```

### Run inference.

The notebooks presented in ```bbh-subpopulations-scripts/inference/``` correspond to 7 of the 8 population inference analyses presented in this paper. To reproduce the analyses (instructions for the 8th at the bottom):

#### GWTC-3, (Analysis corresponding to Figure 1 of the paper)

```
$ mkdir -p /path/to/run-dir-gwtc3
$ cd /path/to/run-dir-gwtc3
$ cp /path/to/bbh-subpopulations-scripts/inference/gwtc3-clean.ipynb .
$ cp /path/to/bbh-subpopulations-scripts/inference/gwtc-3.yaml .
```
Edit the ```.yaml``` file by replacing all occurences of ```/home``` to ```/path/tp/data-dir/home```

Run each cell in the notebook, untill you reach the section labeled "Run-inference".

Below that section is a block of python code. Paste that into a script named ```run_gppop.py``` in the same directory and then:

```
$ python run_gppop.py
```
This code will generate a netcdf file named ```gppop_posterior_*.nc``` which contains the hyper-posterior samples. Optional, run the remaining cells of the notebook for GP-related sanity checks.



#### Simulation study, (Analysis corresponding to Figures 2 and 3 of the paper)

Example given for ```sims_uncorrelated.ipynb```, exactly same for other notebooks:

```
$ mkdir -p /path/to/run-dir-sims_uncorrelated
$ cd /path/to/run-dir-sims_uncorrelated
$ cp /path/to/bbh-subpopulations-scripts/inference/sims_uncorrelated.ipynb .
$ cp /path/to/bbh-subpopulations-scripts/inference/plpp_mdc.yaml
```
Edit the ```.yaml``` file by replacing ```/home``` to ```/path/tp/data-dir/home```

Edit the following paths in the notebook:

In the cell below the heading "Prepare-Input", change ```/home``` to ```/path/tp/data-dir/home```

Find the cell containing ```config_file = '/home/anarya.ray/gppop-prod/pairing/plpp_mdc.yaml'``` and change that to ```config_file = 'plpp_mdc.yaml'```.

In the next to next cell, again change ```/home``` to ```/path/tp/data-dir/home```

Run the rest of the cells untill you reach the section labeled "Run-inference".

Below that section is a block of python code. Paste that into a script named ```run_gppop.py``` in the same directory and then:

```
$ python run_gppop.py
```
This code will generate a netcdf file named ```gppop_posterior_*.nc``` which contains the hyper-posterior samples. Optional, run the remaining cells of the notebook for GP-related sanity checks.

**The instructions for generating all the other ```sims_*``` notebooks are exactly the same, with changed names **


#### GWTC-3, appendix A, (Analysis corresponding to Figures 4,5,6 of the paper)

Reproduce the analysis of Figure 1 mentioned above twice, by first all occurences of ```kappa=3.0``` to ```kappa=1.1``` and then by changint to ```kappa=4.4```.

**Note** Don't forget to change the names of the netcdf file in the python script to avoid overwriting.

#### GWTC-3, appendix B, (Analysis corresponding to Figures 7,8,9 of the paper)


```
$ mkdir -p /path/to/run-dir-gwtc3
$ cd /path/to/run-dir-gwtc3
$ cp /path/to/bbh-subpopulations-scripts/inference/gwtc3-clean_vary_bins.ipynb .
$ cp /path/to/bbh-subpopulations-scripts/inference/gwtc-3.yaml .
```
Edit the ```.yaml``` file by replacing all occurences of ```/home``` to ```/path/tp/data-dir/home```

Run each cell in the notebook, untill you reach the section labeled "Run-inference".

Below that section is a block of python code. Paste that into a script named ```run_gppop.py``` in the same directory and then:

```
$ python run_gppop.py
```
This code will generate a netcdf file named ```gppop_posterior_*.nc``` which contains the hyper-posterior samples. Optional, run the remaining cells of the notebook for GP-related sanity checks.


### Reproducing the plots
```
$ mkdir -p /path/to/plot-dir
$ cd /path/to/plot-dir
$ cp /path/to/bbh-subpopulations-scripts/curate_plots.ipynb .
```

Below the headings with names of each analysis, set ```run_dir=``` with the path to the directory in which the analysis was run/ in which the ```.nc``` file was generated.

Then execute all cells of the notebook to generate the plots.

### Optional: Reproducing simulated datasets