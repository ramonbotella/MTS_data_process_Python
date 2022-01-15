#! python3
# cyclictestMTSuniversal - A Python script to process
#                          data from cyclic testing
#                          from an MTS multiuprupose loading frame.
#
###########################################################################
##
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
##        data file structures and returns two result files with user-choice
##        name '.csv’ and 'xlsx'. In addition the script plots the fatigue curve
##        and finds the failure point by detecting the maximun curvature
##        point on the second part of the curve using kneedle package. The '.xlsx'
##        contains the main parameters of interest, the plot of the fatigue curve
##        and the failure point coordinates.
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
import matplotlib.pyplot as plt
import os
from kneed import KneeLocator


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

while True:
        try:
                side1 = float(input('Enter dimension 1 in mm: '))
        except ValueError:
                print("""Sorry, seems like you didn't type a number. Please, try again:""")
                print('###############################################################')
                continue
        else:
                break

while True:
        try:
                side2 = float(input('Enter dimension 2 in mm: '))
        except ValueError:
                print("""Sorry, seems like you didn't type a number. Please, try again:""")
                print('###############################################################')
                continue
        else:
                break

crossSection = side1*side2
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

# Create dir to store png of hysteresis loops
##os.mkdir('loopAnimation')

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
        stressamp = forceamp*9.81/crossSection
        stressmin = forcemin*9.81/crossSection
        stressmax = forcemax*9.81/crossSection

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
        cycleData = cycleData + [[datacycle[0,(numberRows-1)],
                                  dispamp, forceamp,stressamp,
                                  loop,phase,dispmax,dispmin,
                                  forcemax,stressmax,forcemin,
                                  stressmin,60*(dispmax-dispmin)/0.05]]
##        plt.plot(strain,force)
##        plt.xlabel('Displacement (mm)')
##        plt.ylabel('Force (kgf)')
##        plt.title('Evolution of the hysteresis loop with the number of cycles')
##        plt.savefig('loopAnimation/loop'+str(int(datacycle[0,(numberRows-1)]))+'.png')
##        plt.clf()
# Stores data into a specimenprocess.dat    
cycleData = np.asarray(cycleData)                                  
               
# Create a list with name of columns of results file
colNames = ['Cycle','Disp amp.(mm)','Force amp. (kgf)','Stress amp. (MPa)',
            'Loop area (mm*kgf)','Phase (º)','Disp. max. (mm)',
            'Disp. min. (mm)','Force max. (kgf)','Stress max. (MPa)',
            'Force min. (kgf)','Stress min. (MPa)','Loading rate (mm/min)']

# Convert results to a pandas dataframe and store it in a csv file
df_res = pd.DataFrame(data=cycleData, columns=colNames) 

# Record the cross-section area in the last column, first row
df_res['Cross-section area (mm^2)'] = ''
df_res.iloc[0,-1] = crossSection

# Extracts the 2nd half of the data to calculate failure point
# as the point of max curvature

X = df_res.iloc[int(df_res.shape[0]/2):,0]
Y = df_res.loc[X.index[0]:,'Disp. min. (mm)']

# Max curvature locator
kneedle = KneeLocator(X,Y, S = 1.0, curve = 'concave', direction = 'decreasing')

# Find Y for knee C point
failureCycle = round(kneedle.knee)
f = np.polyfit(X,Y,9)
p = np.poly1d(f)
failureDisp = p(failureCycle)

# Stores the key parameters of the test in first row and last columns of the dataframe
df_res['Failure cycle'] = ''
df_res.iloc[0,-1] = failureCycle

df_res['Mean Force max-min (kgf)'] = ''
df_res.iloc[0,-1] = 2*df_res['Force amp. (kgf)'].mean()

df_res['Mean loading rate (mm/min)'] = ''
df_res.iloc[0,-1] = df_res['Loading rate (mm/min)'].mean()

df_res['Mean Stress max-min (MPa)'] = ''
df_res.iloc[0,-1] = 2*df_res['Stress amp. (MPa)'].mean()

# Input specimen reference to name results files
fileNameRes = str(input('Enter specimen reference:'))

# Stores dataframe in csv file
df_res.to_csv(fileNameRes+'.csv',index = False, sep=';')

# Stores dataframe in xlsx file, plots fatigue curve and failure point

writer = pd.ExcelWriter(fileNameRes+'.xlsx', engine='xlsxwriter')

df_res.to_excel(writer, sheet_name='Processed_data', index=False)
workbook = writer.book
worksheet = writer.sheets['Processed_data']

worksheet.write('S1','Failure Displacement (mm)')
worksheet.write('S2',failureDisp)

chart = workbook.add_chart({'type':'scatter'})

chart.add_series({'categories':'=Processed_data!A2:A10000',
                  'values':'=Processed_data!H2:H10000',
                  'name': fileNameRes
                   })

chart.add_series({'categories':'Processed_data!$O$2',
                  'values':'=Processed_data!$S$2',
                  'name': 'Failure',
                  'marker': {'type': 'square',
                             'size': 8,
                             'border': {'color': 'black'},
                             'fill':   {'color': 'red'}},
                  'data_labels': {'value': False,
                                  'category': True,
                                  'num_format': '#,##0',
                                  'border': {'color': 'red'},
                                  'fill':   {'color': 'yellow'}},
                   })
chart.set_x_axis({
    'name': 'Cycle',
    'name_font': {'size': 14, 'bold': False},
    'num_font':  {'italic': False },
    'num_format': '#,##0'
})

chart.set_y_axis({
    'name': 'Min. Displacement (mm)',
    'name_font': {'size': 14, 'bold': False},
    'num_font':  {'italic': False },
    'crossing':'min',
    'num_format': '#,##0.0'
})

chart.set_size({'width': 720, 'height': 400})

worksheet.insert_chart('C10',chart)

workbook.close()


print("""The results have been stored in """+ fileNameRes +""".csv and """+ fileNameRes +""".xlsx files""")
print("The failure cycle is: " + str(failureCycle))

# Plots the shape of the 2nd part of the curve and the failure point
plt.plot(X,Y,'-',failureCycle,failureDisp,'+')
plt.title('Curvature analysis to find failure cycle')
ax = plt.gca()
ax.set_facecolor((0.898, 0.898, 0.898))
fig = plt.gcf()
plt.xlabel('Cycles')
plt.ylabel('Displacement (mm)')
plt.text(failureCycle, failureDisp, str(failureCycle))
plt.show()


df_res.plot(x='Cycle',y='Disp. min. (mm)')
plt.show()
