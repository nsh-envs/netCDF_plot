# This repository is used for plotting and analyzing netCDF files in Python.
This repository contains some usefull notebooks for working with netCDF files. It also showcases easy ways to plot with various kinds of selection types.

Original author of the notebook labeled 'zhe' is Zhuyun Ye, Postdoc at ENVS, Aarhus Univeristy.

## Files and directories.
### framesX
These directories contain frames for making animations with the DEHM output data. The number *X* refers to the domian of the model.

### make\_video
Small bash script for compiling a video from a directory of frames.

### netcdfX
These directories contains the data. However it can be large amounts. Thus, they will (at some point) contain links to a frozen archive at my university's online filesystem *ERDA*. *TODO: The notebook should also contain a call to wget for automatically downloading the files based on the 'files.txt' file in each directory.*

### plot\_zye.ipynb
This file contains the original file i got from the original author.

### plot\_netCDF.ipynb
This includes my own modification and example code for plotting various netCDF datasets as well as creating animations for DEHM output.

### video/
This is where animations will be stored.
