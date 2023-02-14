# parameters important to the peak finding and ID can be easily changed in this file

common_isotopes = ["Cs-137", "Co-60", "Cs-134", "Bi-207", "Bi-214", "Pb-207", "Pb-210", "Pb-214", "Sr-85", "Ag-110m"]
#common_isotopes = ["Ir-192"]

# for standard
common_isotopes = ["Ba-133", "Co-57", "Ce-139", "Sr-85", "Cs-137", "Mn-64", "Y-88", "Zn-65"]


# detector/spectrum specific:
#c, M = -12.236, 0.9085    # linear energy calibration parameters

# gamma spectrometer fk
c, M = 0.113, 0.1794

#a, b, c = 4.9570877849999996e-05, 0.8177533084999999, 11.660699480000002    # parabolic energy calibration parameters
FWHM = 5   # estimate for peak full width at half maximum
# scintillator = "LaBr3(Ce)"    # scintillator material
# size = 20     # scintillator thickness

scintillator = "HPGe"    # scintillator material
size = 61    # scintillator thickness
time = 1000  # time over which spectrum accumulated

# parameters to adjust sensitivity of peak finding, lowering confidence and intensity threshold will find
# peaks more easily but may identify some non-peak features as peaks additionally:
confidence = 3  # confidence factor, f, for second difference method of peak finding
intensity_threshold = 1E-3 # old
intensity_threshold = 0.001  # threshold for height of detected peaks in %# threshold for height of detected peaks in %
LLD = 60  # lower channel number threshold, everything below will be ignored
HLD = 8192  # higher channel number threshold, everything above will be ignored

# isotope ID parameters:
energy_window_size = 5  # number of keV to include isotopes with energies +/- of detected peak energy
category="All"