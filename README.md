# simGRO
Monte Carlo simulation of RNA polymerase II pausing and elongation

## Quickstart
- Disclaimer: this only works in Mac OS X. Kept getting seg faults in ubuntu.

### Download and install
cd to your installation directory.
```
git clone
cd simGRO
make
```
Ignore all warning messages.
Makes two executable files.
- simGRO : command line GRO-seq simulator
- GROgu : GRO-seq simulator graphical UIed.

### simGRO
Run simGRO. You will need all the parmeter files in the 'default' directory.
```
./simGRO -p default/pref.txt -o result/default.result.txt
```
Run plotGRO.R script to generate the plot of the simulated result.
```
Rscript plotGRO.R
````
![alt txt](https://github.com/h-kwak/simGRO/blob/main/result/default.plot.PNG)

### GROgu
GUI version of simGRO. Only works in Mac OS X.
Start running in terminal commandline. Need the 'gu' directory
```
./GROgu
```
![alt txt](https://github.com/h-kwak/simGRO/blob/main/gu/gu.capture.PNG)

Right click to activate pop-up menu. Keyboard shotcut keys also available (ESC to quit).

## Introduction
Transcription is a dynamic process. On and off rates of transcription factor (TF) binding, RNA polymerase II (Pol II) pausing and elongation, nascent RNA cleavage and termination shape the landscape of Pol II traversing on a gene and dictate the amount of RNA production. Through simulating Pol II dynamics, the following questions will be addressed.
1. Determine the effect of each step, in particular the pausing escape step.
2. Devise experimental parameters that accurately reflect each step of transcription.
3. Train a deep learning machine that determines which step is regulated in a gene.

## SimGRO

### Usage
```
Arguments
-p      :       Preference file name
-o      :       Output file name
-e      :       Equillibrium time (default=600s)
-s      :       Equillibrium speed (default=5x)
-b      :       Bin size (default=10)
-r      :       Simulation resolution (default=50 ms)
```
- Preference file: file containing gene parameters and kinetic constants. Further breakdown on the next section.
- Output file: tab delimited file of recording time, start position, end position, read count per DNA copy in long format dataframe. Use the 'plotGRO.R' R script to produce plot a profile.
- Equillibrium time: time needed to reach equillibrium of Pol II on the gene. Generally, 4x longest half life of the bindings constants or genelength divided by elongation rate, whichever is greater.
- Equillibirum speed: factor of equillibration speed relative to the main simulation. Equillibration can approximated and ran faster.
- Bin size: size of read count bins in the output file in base pairs.  
- Simulation resolution: delta T per each simulation cycle. Note that simGRO uses pre-calculated poisson time function for transcription events, and can handle relatively long delta T. This speeds up the simulation while maintaining accuracy.

### Preference file breakdown
Default preference file at `default/pref.txt`

```
# Default model
```
Lines starting with hash(#) are for annotation purposes

```
# Gene length and positions
GeneLength	10000
PostPolyALength	1800
PromoterLength	200
PauseSite	40
```
Gene lengths and positions. Labels are case insensitive, followed by a tab and the number in integers. Promoter starts at 0 bo, TSS is at 200 bp, pause site at 240 bp, and gene is 10 kb long, followed by 1800 bp flanking region.

```
# Initiation rate limiting transcription factor occupancy parameters (once every 500s, stays for 20s)
# Rows are the function of time
TF_on	default/tf_on.txt
TF_off	default/tf_off.txt
```
Parameter files defining transcription factor binding and dissociation rates. File name(default/tf_on.txt) follows the labels after a tab. Files are tab delimited table of 2 variable function, usually Pol II position on the x (column) and time on the y (row). TF binding is only time dependent, and has 1 column. Further parameter file breakdown on the next section. TF binding is the first rate limiting step, set to happen once every 500 seconds. TF binding is relatively dynamic, dissociating after the average of 20 seconds. Default values modeled after TATA Box binding protein (TBP).

```
# Pol II recruitment rate, rows are the function of time
Recruitment	default/recruitment.txt
```
Pol II recruitment rate function. This is the non-rate limited intrinsic Pol II recruitment rate.

```
# Pausing factor parameters, function of time 
PF_on	default/pf_on.txt
PF_off	default/pf_off.txt
```
Pausing factor recruitment and escape rates. Once Pol II reaches the pausing site, fast PF_on rate captures Pol II. PF_off is the second rate limiting step.

```
# Elongation rate, function of position and activity
Elongation	default/elongation.txt
```
Elongation speed set at 30 bp/sec constant. Elongation speed can be defined at different rates depending on positions (columns). Pol II 'activity' is an intrinsic property of the activity of Pol II molecule (such as Ser2 phosphorylation, on rows) that can affect the speed.

```
# RNA cleavage rate (polyA), by position and activity
Cleavage	default/cleavage.txt
```
Cleavage polyadenylation rate, dependent on the position. Is 0 on the gene body up to 10 kb (+200 bp), and goes up after the PAS. Can further be set to be dependent on Pol II activity.

```
# Post-cleavage elongation rate, function of position and activity
PostPA	default/postpa.txt
```
Post-poly(A) elongation rate, set at 10 bp/sec. Post-poly(A) rate applies only after the RNA cleavage event. 

```
# Termination rate, function of position and activity
Termination	default/term.txt
```
Termination rate. Note that the default cleavage, postPA elongation, and termination rates are arbitrary values to recapitulate typical post-PAS transcription profiles and not affect gene body profiles.

```
# Pol II initial activity distribution (scale of 1 to 100, function of time)
Activity	default/activity.txt
```
Pol II activity definition. Currently not used, and all Pol II have 100% activity. Reserved to model activity dependent Pol II behaviors.

```
# Recording timepoints (seconds)
Record	default/record.txt
```
File containing simulated profile recording timepoints. Default recording at 0 sec and 60 sec.

```
# DNA copies to simulate
Copies	2000
```
Number of DNA template to simulate. Can be adjusted to your system performance and simulation depth.

```
# End of preference definition
```

### Parameter file breakdown
TF`default/tf_on.txt`
```
TFon_R)time	 0
0            0.002
600          0.002
```
