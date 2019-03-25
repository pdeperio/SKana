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
```
cd ana
root plot.C
```
