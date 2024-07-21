# bbh-subpopulations-scripts

Scripts for running the binned gaussian process analyses presented in [2404.03166](https://arxiv.org/abs/2404.03166). To run these scripts the following packages are necessary:

- https://github.com/AnaryaRay1/bbh-subpopulations-scripts/tree/main ([9f853ce](https://github.com/AnaryaRay1/gppop/commit/9f853cecc3ab1c5ee4358baaf6876492d59c1945))
- https://git.ligo.org/reed.essick/gw-distributions.git ([5d3804](https://git.ligo.org/reed.essick/gw-distributions/-/commit/5d38045f93cf5367ff031689ea27e1d640f0f3c1))
- https://pypi.org/project/bilby/ (latest)
- https://pypi.org/project/gwpopulation/ (latest)

as well as the datasets uploaded to: https://zenodo.org/records/12746314

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


