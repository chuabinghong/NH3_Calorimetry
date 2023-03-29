# <u>Calorimetry Analysis</u>
*Repository for analysis code for the calorimetry data of NH3-H2O systems*

All data is saved as CSVs with the first 8 rows being metadata. The plots show the temperature profile, heat flow, heat capacity (even though it is not correct over the entire curve yet due to fractional crystallization), and lastly a Shomate curve fit of the data with peaks removed. 
For the code ('CalorimeterProcessing.py'), you can ask me any questions you like, but this it is pretty basic and was just helpful for me to organize and play with the data. If you just want to look at the data, all you really have to use is 'plots()' explained below.

The 'plot()' functions just allows you to plot all the samples or just one. You either pass 'all' (default) or the number of the wt.% sample you're interested in. (I.e., plots(4) would return the plot for the 4 wt.% data). You will find the call to 'plot()' at the end of the python file. 
This is really just to prevent a bunch of plots printing every time if its not needed.

We noticed that the 'double' peaks somewhat intensified with increased wt.% and did not show up at all in the 2 wt. % sample. The shape of the double peaks, however, was not consistent throughout the runs.

Jack

*For future contact please use <b>jdiab@ucla.edu</b>*
