# Bayesian adaptive dual control of deep brain stimulation in a computational model of Parkinson's disease

##### Logan L Grado, Matthew D Johnson, Theoden I Netoff

The code to run the basal ganglia thalamocortical system (BGTCS) mean field model (MFM) is provided here. Additionally, we have provided the simulated continuous deep brain stimulation (cDBS) and phasic DBS (pDBS).

## Installation
1. Create and activate virtual environment
```
$ virtualenv pyenv
$ source pyenv/bin/activate    # OSX/Linux
$ pyenv\Scripts\activate.bat   # Windows
```
> **Note:** Using a virtual environment isn't necessary, but will prevent conflicts with any packages you already have installed on your system.


2. Install requirements
```
$ pip install -r requirements.txt
```

## Running the Model
##### The model is contained in `mfm.py`. To run the model, run `mfm.py`.
```shell
$ python mfm.py
Run Info
--------
length   : 50.0 s
DD       : True
cDBS     : False
pDBS     : False

SWIFT Parameters
----------------
f     : 29.0 Hz
tau_s : 0.2397 s
tau_f : 0.0479 s

0:00:05 [==================== 100.00 % ====================] 0:00:00
 
Saving data...
  RunID: 000
```
Each run is automatically numbered and saved in `data/`.

##### For usage instructions and help, use `-h` or `--help`
```shell
$ python mfm.py --help
BGTCS MFM

Usage:
  mfm [options] [<key>=<value>]...
  
Options:
  -l --list    List all options
  -h --help    Show this screen
```

##### To list available options and their defaults, use `-l` or `--list`
```shell
$ python mfm.py --list              
Available options:
Option           Default
---------------  ---------
Cm               0.0001
DD               True
RunID            -1
cDBS             False
cDBS_amp         3.0
cDBS_f           130.0
cDBS_width       60.0
dt               0.001
pDBS             False
pDBS_amp         2.38
pDBS_phase       2.24
pDBS_power_thr   -28.57
pDBS_ref_period  0.3
pDBS_width       60
state_target     p1
stim_start       0.0
stim_target      STN
swift_c          10
swift_f          29
swift_s2f        5
swift_tau_s      0.2397
tstop            50.0
verbose          True

For more details, such as units, etc, look at the source of MFM._load_params().
```

##### Examples:
```shell
$ python mfm.py DD=False                # turn DD off
$ python mfm.py tstop=100 cDBS=True     # run for 100s, turn cDBS on
$ python mfm.py pDBS=True               # turn pDBS on
```

## Plotting the Results
By default, all model runs are saved in `data/` as `<RunID>.mfm`. The runs are saved simply by pickling the MFM object. The script `plot.py` is provided to re-load and plot a run after saving and exiting. 

##### Plot run
```shell
$ python plot.py 0
Run Info
--------
length   : 50.0 s
DD       : True
cDBS     : False
pDBS     : False

SWIFT Parameters
----------------
f     : 29.0 Hz
tau_s : 0.2397 s
tau_f : 0.0479 s
```

`plot.py` either takes a RunID `int` or filename `str`. It also prints metadata about the run.

## Figure 3

To generate figure 3 from the paper, run `fig_3.py`.
```shell
$ python fig_3.py
```
