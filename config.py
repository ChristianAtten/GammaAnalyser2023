# parameters important to the peak finding and ID can be easily changed in this file

common_isotopes = ["Cs-137", "Co-60", "Cs-134", "Bi-207", "Bi-214", "Pb-207", "Pb-210", "Pb-214", "Sr-85", "Ag-110m"]

# detector/spectrum specific:
c, M = 0, 1     # energy calibration parameters
FWHM = 14    # estimate for peak full width at half maximum

# changed here
scintillator = "LaBr3(Ce)"    # scintillator material
size = 20     # scintillator thickness
time = 300  # time over which spectrum accumulated

# parameters to adjust sensitivity of peak finding, lowering confidence and intensity threshold will find
# peaks more easily but may identify some non-peak features as peaks additionally:
confidence = 3  # confidence factor, f, for second difference method of peak finding
intensity_threshold = 0.001  # threshold for height of detected peaks in %
LLD = 0  # lower channel number threshold, everything below will be ignored
HLD = 1e9  # higher channel number threshold, everything above will be ignored

# isotope ID parameters:
energy_window_size = 5  # number of keV to include isotopes with energies +/- of detected peak energy
category="All"