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
RNA polymerase II (Pol II) is a dynamic enzyme carrying out the key step in transcription. Transcription factors (TFs) bound to promoters recruit Pol II, but Pol II pauses soon after initiating RNA synthesis. Paused Pol II escapes to productive elongation making nascent RNA molecules, until the RNA is cleaved and Pol II is terminated. The on and off rates of TFs and pausing factors, as well as Pol II elongation speed and termination rates shape the overall Pol II landscape on a gene. How these steps dictate the activity of RNA synthesis has been the major focus of the studies of transcription mechanisms.

