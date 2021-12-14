# MTS_data_process_Python
The software from the MTS hydraulic press employed
for asphalt materials testing at the Road Reseach Laboratory
of the UPC-BarcelonaTech acquires detailed data every certain
number of points (typically, every 100 cycles). For each
cycle acquired, the software registers approximately 50 data
points (1 data point each (freq.*50)^1 seconds). This produces
large data files where meaningful data is partitioned in blocks
of 50 data points, with blank and descriptive rows containing
strings in between. The fact that the number of data points
acquired in each cycle fluctuates between 49 and 51,
complicates the automatization of computing the main parameters
that describe each cycle (force, displacement amplitudes and
delay between the two of them, i.e., phase angle).

This script automatizes this process. This script is tailored
to the most common data file structures in the Road Reseach
Laboratory  of the UPC-BarcelonaTech, but it can be easy adapted
to different column distributions, names, units, etc,…

In the current version this script can deal with two different
data file structures and returns the ‘specimenprocess.dat’ file.
The result file is a comma separated file where the following
values for each recorded cycle are reported: Cycle number,
piston displacement amplitude in mm, force amplitude in kgf,
hysteresis loop area in mm·kgf, phase angle in the first
quadrant in degrees , maximum and minimum displacement during
the cycle in mm, maximum and minimum force registered during
the cycle in kgf, and loading rate in mm/min (speed of the
piston peak-valley).

Sample data files with the two data structures the script
accepts can be found on the following repository:
https://github.com/ramonbotella/MTS_data_process_Python/
