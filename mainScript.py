import pandas as pd
import sys
import time
import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import analyzingFunctions
import spectrometryConfigurations
import plottingfunctions
from datetime import datetime
import scipy

########################################################################################################################
# 1. Define Path and File of Spectrum
########################################################################################################################

#path = 'C:/Users/chris/Dropbox/Universität im Fürstentum Liechtenstein/UFL-Liechtenstein_Datenauswertung/' \
#       'GammaRay-Spectrometry/existingFrameWork_Sam/inputDaten/'
#file = path + 'spectrum.txt'


# path for sam data
#path='C:/Users/chris/Dropbox/Universität im Fürstentum Liechtenstein/UFL-Liechtenstein_Datenauswertung/GammaRay-Spectrometry/GammaAnalyzer/inputData/datenSam/'
#file = path + "spectrum_eu152_2700s.txt"

# path for fk data
path='C:/Users/chris/Dropbox/Universität im Fürstentum Liechtenstein/UFL-Liechtenstein_Datenauswertung/GammaRay-Spectrometry/GammaAnalyzer/inputData/datenFeldkirch/Kalibrationsstandard/'
file = path + "20220401_StandardHiEng_BC9984.csv"
file = path + "highEng_10_10_2022.csv"


path = 'C:/Users/chris/Dropbox/Universität im Fürstentum Liechtenstein/UFL-Liechtenstein_Datenauswertung/GammaRay-Spectrometry/GammaAnalyzer/inputData/datenFeldkirch/Iridium/'
#file = path + "checkcableIridiumWischprobe.csv"



pathData='C:/Users/chris/Dropbox/Universität im Fürstentum Liechtenstein/UFL-Liechtenstein_Datenauswertung/GammaRay-Spectrometry/GammaAnalyzer/'
pathInputData=pathData+'inputData/'



now = str(datetime.today())
todayDate = now.split(' ')[0]
pathOutputData=pathData+'outputData/' +todayDate + '/'
if not os.path.exists(pathOutputData):
    os.makedirs(pathOutputData)

isotopeLib = 'IsotopeLibrary_splitted_IR.csv'
fileIsotope = pathInputData+isotopeLib

# Read spectrum from txt file
data=pd.read_csv(file,delimiter=';',index_col=False)
data_help = data
#data = np.genfromtxt(file, delimiter="\n", dtype=int)   # read spectrum as a numpy array
test=data.iloc[:,[2]]
test2=test.T.to_numpy()
data=test2
data=data[0][1:len(data.transpose())]
#data=np.round(data/1000,0)
x=0
########################################################################################################################
# 2. Set Spectrometry Configurations in spectrometryConfigurations.py
########################################################################################################################
try:
    # from config import FWHM, HLD, LLD, intensity_threshold  # get config values for use in peak finding algorithm
    from spectrometryConfigurations import FWHM, HLD, LLD, \
        intensity_threshold  # get config values for use in peak finding algorithm

except:
    print("FWHM estimate, HLD and LLD must be provided in config.py")
    sys.exit()

# if the discriminators are outside of range of data, can't use them so just use all data
if spectrometryConfigurations.HLD > len(data) or spectrometryConfigurations.LLD < 0:
    print("LLD or HLD not in range of spectrum, discriminators will be disregarded")
    spectrometryConfigurations.HLD = len(data) - 1
    LLD = 0

# empirical values for peak finding algorithm found by Mariscotti
z = 5
w = int(0.6 * spectrometryConfigurations.FWHM)
if (w % 2) == 0:
    w += 1
m = int((w - 1) / 2)

# attempt an efficiency correction if one has been provided
try:
    y = analyzingFunctions.efficiency_correction(pathInputData,data)
    efficiency_corrected = True
except:
    print("Energy-wise efficiency correction missing or invalid, raw spectrum will be used")
    y = data
    efficiency_corrected = False

########################################################################################################################
# 3. Find peak centres and their index in the spectrum
########################################################################################################################
signals = np.zeros(data.shape)
for i in range(FWHM - 20, FWHM + 21):
    # empirical values for peak finding algorithm found by Mariscotti
    z = 5
    w = int(0.6 * i)
    if (w % 2) == 0:
        w += 1
    m = int((w - 1) / 2)
    if i < 4:
        continue
    signals += analyzingFunctions.peak_finding(y, z, m, i)
locs = np.flatnonzero(signals).tolist()

########################################################################################################################
# 4. Remove linear background
########################################################################################################################
x, removed = analyzingFunctions.remove_bg(y, locs, FWHM)

########################################################################################################################
# 5. Remove detected peaks which are relatively small or outside of discriminators
########################################################################################################################
locs2 = locs.copy()
for i in range(len(locs)):
    if locs[i] <= LLD or locs[i] >= HLD:
        locs2.remove(locs[i])
        continue
    if y[locs[i]] / max(y[LLD:HLD]) < intensity_threshold * 100 or locs[i] < 2 * FWHM + 1:
        locs2.remove(locs[i])
        continue
    for j in range(len(locs2)):
        if locs[i] in range(locs2[j] - FWHM, locs2[j] + FWHM) and locs[i] != locs2[j]:
            locs2.remove(locs[i])
            break

#locs2=[238, 421,1138, 1650, 1720, 1765, 2609,3412]
# fit the peaks with gaussian and return goodness of fits, areas and FWHMs
new_locs, goodness_of_fits, fits, net_areas, gross_areas, FWHMs = analyzingFunctions.peak_fitting(x, locs2, FWHM, removed)

# all located peak centre channels converted to energies
try:
    energies = analyzingFunctions.linear_channel_to_energy(new_locs)

    # new
    testX = data_help.iloc[new_locs,1]
    testX2 = testX.T.to_numpy()
    energies = testX2

except ImportError:
    energies = analyzingFunctions.parabolic_channel_to_energy(new_locs)

# attempt to identify sources responsible for spectrum observed
table_of_sources = analyzingFunctions.source_lookup(fileIsotope,energies, goodness_of_fits, gross_areas, efficiency_corrected)

# all obtained data to a dataframe in order to be plot in a table then returned with the fits for all peaks
dframe = pd.DataFrame(table_of_sources, columns=["Isotope", "Energy (keV)", "RI (%)", r"$\bar{R}^2$ of fit",
                                                     "CPS", "ACF"])

########################################################################################################################
# 6. Plot outcome
########################################################################################################################
plottingfunctions.plotOutcome(pathOutputData, file,data, dframe, fits, FWHMs, net_areas, gross_areas)