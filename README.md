# simGRO
Monte Carlo simulation of RNA polymerase II pausing and elongation

## Quickstart
- Disclaimer: this only works in Mac OS X. Kept getting seg faults in ubuntu. Can't figure out why.

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
Start running in terminal commandline. Need the 'gu' directoru
```
./GROgu
```
