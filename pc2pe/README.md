# pc2pe: Relative gain analysis

## Usage:

1. Checkout repository:
```
git clone https://github.com/pdeperio/SKana.git
```

2. Compile on sukap:
```
source SKana/env.sh
cd SKana/pc2pe
make
```

3. Run:

* Specify datasets and input directory coming from ```tqreal``` processing in ```sub_job.sh```
* Submit jobs with: ```./sub_job.sh```

4. Concatenate histograms:
See ```hadd.sh``` for example.

5. Analyze:

* Scripts can be found in:
```
cd ana
```

* Then create pc2pe TTree, including diagnostic plots, like hit-time window selection:
```
root plot.C
```

* Double check time stability with:
```
root plot_time.C
```
(If you want to see entire time distribution of all runs, need to modify runs processed in `sub_job.sh` and `hadd.sh` above.)

* Make grouping plots (manually select datasets to overlay within):
```
root plot_position.C
```

* And pc2pe 1D (2D) projections (correlations):
```
root plot_pc2pe.C
```

* Append SPE analysis information to TTree (warning, this overwrites the file produced in `plot.C`): 
```
root -b -q 'plot_spe.C(1)'
```

* Then make SPE plots and correlations:
```
root plot_spe.C
plot_spe_peaks.C
```

* Finally, once everything looks good, prepare the relative gain correction tables with:
```
root make_pc2pe.C
```
