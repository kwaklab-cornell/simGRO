# simGRO
Monte Carlo simulation of RNA polymerase II pausing and elongation

## Quickstart
- Disclaimer: this only works in Mac OS X. Kept getting seg faults in ubuntu.

### Download and install
cd to your installation directory.
```
git clone https://github.com/h-kwak/simGRO.git
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
Lines starting with the hash (#) are for annotation purposes

```
# Gene length and positions
GeneLength	10000
PostPolyALength	1800
PromoterLength	200
PauseSite	40
```
Gene lengths and positions. Labels are case insensitive, followed by a tab and the number in integers. Promoter starts at 0 bp, TSS is at 200 bp, pause site at 240 bp, and gene is 10 kb long, followed by 1800 bp flanking region.

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
TF binding parameters at `default/tf_on.txt`. (Note that the tabs are expanded)
```
TFon_R)time	 0
0            0.002
600          0.002
```
First line is a header of x paramters. TF binding is not dependent on Pol II position (x), and only one value of 0.
Rows are the y parameter values (time), followed by the TF binding rate constant (per second). Rate constant of 0.002 is once every 500 seconds.
This was set to be independent of time, but can be adjusted to model dynamic changes.
- TF_off, PF_on, PF_off parameter files follow the same rule.

Elongation rate parameters at `default/elongation.txt`.
```
ElongationRate_C)distance:R)activity  0    10000
0                                     30   30	
100                                   30   30
```
Elongation rate is a two variable function (Pol II position on the x column, and Pol II activity on the y row). Default is set at 30 bp/second, but can be varied. Values are interpolated and extrapolated linearly.

Cleavage rate parameters at `default/cleavage.txt`.
```
CleavageRate_C)pos:R)activity  0   10199   10200   12000
0                              0   0       0.05    0.05
100                            0   0       0.05    0.05
```
Cleavage rate is also defined similarly. Note that there is an abrupt increase from 0 between 10199 and 10200. Cleavage and termination constants are set up arbitrarily to match a typical post-PAS Pol II distribution.

### plotGRO.R
Change the following part to match your result file, and run `Rscript plotGRO.R` to generate the plots.
```
# Simulation information
geneLength = 10000
postPolyALength	= 1800
promoterLength = 200
pauseSite = 40
resultFile = "result/default.result.txt"
outFile = "result/default.result.pdf"
```

## GROgu
GROgu is a GUI version of simGRO. The difference from simGRO is that GROgu use constant rate, easier to edit the parameter files for simpler illustrative purposes. 

### Usage
```
./GROgu <preference file(optional)>
```
Without specifying the preference file, GROgu will use `gu/pref.txt`.
The graphical simulation will start. Right click to start a pop-up menu. Details on the commands will be updated.
Usefull hotkeys are
- ESC: quit
- Number keys: shift to conditions set on preference file.
Also shift between different conditions through the pop-up menu.

### preference file
```
7
hsp70 NHS
hsp70 HS
Class 2 low esc
Class 2 high esc
Class 1
Class 4
Fast TBP dynamics
gu/par.hsp70NHS.txt
gu/par.hsp70HS.txt
gu/par.class2LowEsc.txt
gu/par.class2HighEsc.txt
gu/par.class1.txt
gu/par.class4.txt
gu/par.fastTBP.txt
```
First line is the number of parameter sets (7 in this case). Each set represents a condition.
Next 7 are the description texts of the conditions. Then the list of parameter files follow. You can assign any other file as long as it follows the parameter file format in the next section.

### parameter files
A parameter file simulating hsp70 gene under non-heatshock condition is at `gu/par.hsp70NHS.txt`.
```
# Hsp70-like gene BEFORE heat shock

# Initiation complex occupancy parameters (once every 500s, stays for 20s)
InitComplex_on	0.002
InitComplex_off	0.05

# Pol II initiation rate
Recruitment	1

# Pause escape rate (once every 500s)
PauseFactor_off	0.002

PostPA-Elongation	15
Termination	0.2
```
The paramters are constant values in GROgu, rather than multivariable functions in simGRO.


## Exercise 1
Regardless of the pause occupancy level, Pol II pausing regulate productive elongation by limiting transcriptional bursting.

Transcription intiation was long assumed to be the rate limiting step of RNA synthesis. Pol II pausing recently emerged as another rate limiting step. Determining which is the limiting step of a gene becomes important, since each involves different factors. We will simulate the condisions where changing the pause escape from low to high (and other parameters remain the same) result in increased gene body Pol II densities.  

### Snapshot
Class 2 genes (active and paused genes) in GROgu default pref.
- Low escape from pause
![alt txt](https://github.com/h-kwak/simGRO/blob/main/gu/gu.c2le.PNG)
- High escape from pause
![alt txt](https://github.com/h-kwak/simGRO/blob/main/gu/gu.c2he.PNG)

Gene body density is apparently higher in high escape, while other parameters are the same and pause site occupancy in low escape is ~6%.

### Simulation
Run R console (or Rstudio), and load
```
source("beskar.R")
```

Set limits for the TF_on, TF_off, and PIC recruitment paramaters
```
par.lim = matrix(log10(c(0.0002, 1,              # TF_on limits: 1 sec to ~1 hr
                         0.005, 0.5,             # TF_off limits: 2 sec to ~3 min
                         0.01, 1)),              # PIC recruitment rate: 1 sec to ~1.5 min
                 ncol = 3,
                 dimnames = list(c("min", "max"),
                                 c("tf_on", "tf_off", "rec")))
```

Shotgun strategy to cover the parameter limits by iterative random sampling. 
```
par = NULL
while(true) {
    par = jetpack(par, par.lim)
    save.params(par, "sample/ex1/pars.txt")
}
````
Function jetpack randomly samples the tf_on, tf_off, and rec parameters, and run simGRO under 3 different pause ecape rates (once every 5s, 12.5s, 50s: simGRO preference files are pref.he.txt, pref.me.txt, pref.le.txt). Run until the sampling saturates the parameter space (unlikely) or time permits. Results are saved in the sample/ex1/result directory.

After (or during) the run, retrieve randomized parameters and simulated read densities
```
data = grep.data("sample/ex1/pars.txt")
```
We will measure fold change in gene body density from low to high escape rate change. This is the column data$GBfc.
Generate 3 variable (tf_on, tf_off, recruitment) plots showing GBfc, pause site density in low escape (LE_PP).
```
plot.data(data, var = "GBfc", file = "sample/ex1/GBfc.pdf")
plot.data(data, var = "LE_PP", file = "sample/ex1/PPle.pdf")
```
Gene body fold change
[!alt txt](https://github.com/h-kwak/simGRO/blob/main/result/GBfc.png)
Pause occupancy
[!alt txt](https://github.com/h-kwak/simGRO/blob/main/result/PPle.png)


## Exercise 2
Dynamic modeling of Pol II elongation wave (Jonkers et al, eLife, 2014).

## Exercise 3
Generate training sets of simulated Pol II profiles of the change in initiation rates or the pause-escape rates. Use these sets to evaluate the perfomance of experimental parameters. Develop a deep learning machine that can tell if a gene is regulated at the initiation or the pause-escape steps across cell types. Discover underlying sequence elements and TFs that differentiate intiation regulated vs pause-escape regulated promoter.
