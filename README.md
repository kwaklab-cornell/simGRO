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
Run simGRO
```
./simGRO -p default/pref.txt -o result/default.result.txt
```
Run plotGRO.R script to generate the plot of the simulated result.
```
Rscript plotGRO.R
````
