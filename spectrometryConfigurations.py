import analyzingFunctions


# parameters important to the peak finding and ID can be easily changed in this file

#common_isotopes = ["Cs-137", "Co-60", "Cs-134", "Bi-207", "Bi-214", "Pb-207", "Pb-210", "Pb-214", "Sr-85", "Ag-110m"]
common_isotopes = []
# detector/spectrum specific:
#c, M = -12.236, 0.9085    # linear energy calibration parameters
c, M = 0,1   # linear energy calibration parameters
c, M = -0.102, 0.2441

#a, b, c = 4.9570877849999996e-05, 0.8177533084999999, 11.660699480000002    # parabolic energy calibration parameters
FWHM = 35   # estimate for peak full width at half maximum
scintillator = "LaBr3(Ce)"    # scintillator material
size = 20     # scintillator thickness

### neu

scintillator = "HPGe"    # scintillator material
size = 61     # scintillator thickness
time = 1000  # time over which spectrum accumulated

# parameters to adjust sensitivity of peak finding, lowering confidence and intensity threshold will find
# peaks more easily but may identify some non-peak features as peaks additionally:
confidence = 3 # confidence factor, f, for second difference method of peak finding
intensity_threshold = 1E-9  # threshold for height of detected peaks in %
LLD = 0  # lower channel number threshold, everything below will be ignored
HLD = 8192  # higher channel number threshold, everything above will be ignored

# isotope ID parameters:
energy_window_size = 10  # number of keV to include isotopes with energies +/- of detected peak energy
#category="Homeland Security"

# CA
category="All"


def spectrometryConfiguration(data):
    try:
        from config import FWHM, HLD, LLD, intensity_threshold  # get config values for use in peak finding algorithm
    except:
        print("FWHM estimate, HLD and LLD must be provided in config.py")
        sys.exit()

    # if the discriminators are outside of range of data, can't use them so just use all data
    if HLD > len(data):
        print("HLD not in range of spectrum, will be disregarded")
        HLD = len(data) - 1
    if LLD < 0:
        print("LLD not in range of spectrum, will be disregarded")
        LLD = 0

    # attempt an efficiency correction if one has been provided
    try:
        y = efficiency_correction(data)
        efficiency_corrected = True
    except:
        print("Energy-wise efficiency correction missing or invalid, raw spectrum will be used")
        y = data
        efficiency_corrected = False

    # find peak centres and their index in the spectrum
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

    if all(LLD > elem or elem > HLD for elem in locs) or len(locs) < 1:
        return pd.DataFrame(columns=["Isotope", "Energy (keV)", "RI (%)", r"$\bar{R}^2$ of fit",
                                     "CPS", "ACF"]), np.zeros(y.shape)

