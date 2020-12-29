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

