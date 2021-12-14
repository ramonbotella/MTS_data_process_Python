#! python3
## MTS_cyclic.py - A Python script to process
##                 data from cyclic testing
##                 from an MTS multiuprupose loading frame.
##
###########################################################################
##        Description:
##        The software from the MTS hydraulic press employed
##        for asphalt materials testing at the Road Reseach Laboratory
##        of the UPC-BarcelonaTech acquires detailed data every certain
##        number of points (typically, every 100 cycles). For each
##        cycle acquired, the software registers approximately 50 data
##        points (1 data point each (freq.*50)^1 seconds). This produces
##        large data files where meaningful data is partitioned in blocks
##        of 50 data points, with blank and descriptive rows containing
##        strings in between. The fact that the number of data points
##        acquired in each cycle fluctuates between 49 and 51,
##        complicates the automatization of computing the main parameters
##        that describe each cycle (force, displacement amplitudes and
##        delay between the two of them, i.e., phase angle).
##
##        This script automatizes this process. This script is tailored
##        to the most common data file structures in the Road Reseach
##        Laboratory  of the UPC-BarcelonaTech, but it can be easy adapted
##        to different column distributions, names, units, etc,…
##
##        In the current version this script can deal with two different
##        data file structures and returns the ‘specimenprocess.dat’ file.
##        The result file is a comma separated file where the following
##        values for each recorded cycle are reported: Cycle number,
##        piston displacement amplitude in mm, force amplitude in kgf,
##        hysteresis loop area in mm·kgf, phase angle in the first
##        quadrant in degrees , maximum and minimum displacement during
##        the cycle in mm, maximum and minimum force registered during
##        the cycle in kgf, and loading rate in mm/min (speed of the
##        piston peak-valley).
##
##        Sample data files with the two data structures the script
##        accepts can be found on the repository.
##
###########################################################################

import numpy as np, math as m
import math as m
import pandas as pd
from scipy.optimize import leastsq                                      
from tkinter import filedialog
from tkinter import *

# Opens window dialog to indicate the location and name of the data file
root = Tk()
root.data =  filedialog.askopenfilename(
        title = "Select data file",filetypes = (
                ("dat","*.dat"),("all files","*.*")))                   
root.destroy()                                                         

# User queries regarding format of data file
while True:
        separator = input('Enter the decimal separator used in the datafile (,/.):')
        if separator == ',' or separator =='.':
                break
        else:
                print('!!!!!!!!!!ERROR!!!!!!!!!!')
                print(""" Please, type ',' or '.' """)
                print('#########################')
                continue
while True:
        language = input('Enter language of datafile (e for English/s for Spanish):')
        if language == 'e' or language == 's':
                break
        else:
                print('!!!!!!!!!!ERROR!!!!!!!!!!')
                print(""" Please, type 'e' for English or 's' for Spanish """)
                print('################################################')
                continue

while True:
        numberRows = int(input('Enter number of rows on data file (4/6):'))
        if numberRows == 4 or numberRows == 6:
                break
        else:
                print('!!!!!!!!!!ERROR!!!!!!!!!!')
                print(""" Please, type 4 or 6 """)
                print('#####################')
                continue
while True:
        try:
                frequency = float(input('Enter test frequency: '))
        except ValueError:
                print("""Sorry, seems like you didn't type a number. Please, try again:""")
                print('###############################################################')
                continue
        else:
                break                

omega = frequency*2*np.pi

if language == 'e':
        pointName = 'Points: '
        timeName = 'Time'
        dispName = 'Axial Desplazamiento'
elif language == 's':
        pointName = 'Puntos: '
        timeName = 'Tiempo'
        dispName = 'Channel 1 Displacement'

# Loading data and measuring data file size in rows
df = pd.read_csv(root.data, header = 4, sep='\t',decimal = separator,
                 on_bad_lines='skip', encoding='mbcs')
print('The selected data file has '+str(df.shape[0])+' lines')

# Initialize list to store cycle data
cycleData = []                                                         
                                                           
# Dataframe with the number of data points
# acquired in each cycle
points = df.loc[df[timeName]== pointName]
points = points.astype({dispName: int})

for index,row in points.iterrows():
        # Isolates the data block for corresponding cycle
        datacycle = df.iloc[(index+3):(index+row[1]+3)]
        datacycle = datacycle.replace({separator:'.'}, regex = True)
        datacycle = datacycle.astype(float)
        datacycle = np.asarray(datacycle)                               

        # Comptutation of amplitude, max. and min. values of the cycle
        dispamp = (datacycle[:,1].max()-datacycle[:,1].min())/2         
        forceamp = (datacycle[:,2].max()-datacycle[:,2].min())/2
        dispmin = datacycle[:,1].min()
        dispmax = datacycle[:,1].max()
        forcemax = datacycle[:,2].max()                                  
        forcemin = datacycle[:,2].min()                                 

        # Computation of the hysteresis loop area using the
        # shoe-lace formula aka Gauss determinant formula
        #       Computes closing loop values 
        add = datacycle[-1,1]*datacycle[0,2]                            
        subs = datacycle[0,1]*datacycle[-1,2]

        #       Computes the rest of the terms
        for k in range(datacycle.shape[0]-1):                           
                add = add + datacycle[k,1]*datacycle[k+1,2]
                subs = subs + datacycle[k+1,1]*datacycle[k,2]
        #       Applies the formula to get the loop area        
        loop = 0.5*abs(add-subs)                                        

        # Computation of the phase angle by fitting sinusoidal
        # to force and displacement signals
        #       Extracts time, force and strain as column vectors      
        t = datacycle[:,0]                                              
        force = datacycle[:,2]
        strain = datacycle[:,1]

        #       Sets initial guesses for the fitting parameters
        force_guess_mean = np.mean(force)                               
        force_guess_amp = forceamp
        force_guess_phase = 0
        strain_guess_mean = np.mean(strain)
        strain_guess_amp = dispamp
        strain_guess_phase = 0

        #       Defines the sinusoidal functions to optimize
        #       by least squares method and fit the parameters
        #       Force signal
        optimize_func_force = (
                lambda x: x[0]*np.sin(omega*t+x[1])+x[2]-force)         
        force_est_amp, force_est_phase, force_est_mean = leastsq(
                optimize_func_force,
                [force_guess_amp,
                 force_guess_phase,force_guess_mean])[0]                                                           

        #       Displacement signal
        optimize_func_strain = (
                lambda x: x[0]*np.sin(omega*t + x[1])+x[2]-strain)                             
        strain_est_amp, strain_est_phase, strain_est_mean =leastsq(
                optimize_func_strain,
                [strain_guess_amp,
                 strain_guess_phase,strain_guess_mean])[0]              

        #       Computes the phase lag between force and displacement
        phase = abs(force_est_phase - strain_est_phase)*360/(2*np.pi)   
        if 90 < phase < 180:
                phase = 180 - phase

        # Appends the data from the cycle to the cycleData list        
        cycleData = cycleData + [[                                       
                datacycle[0,(numberRows-1)],dispamp, forceamp,          
                loop,phase,dispmax,dispmin,forcemax,forcemin,
                60*(dispmax-dispmin)/0.05]]
                                                                        
# Stores data into a specimenprocess.dat    
cycleData = np.asarray(cycleData)                                  
np.savetxt('specimenprocessed.dat',cycleData, delimiter=',',newline='\r\n',
           header="""Cycle,Disp amp.(mm),Force amp. (kgf),\
           Loop area (mm·kgf),Phase (º),Disp. max. (mm),Disp. min. (mm),\
           Force max. (kgf),Force min. (kgf),Loading rate (mm/min)""")                      

# Prints a remainder of where the data has been stored
print("""The processed data has been saved into the file 'specimenprocessed.dat'""")

