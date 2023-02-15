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
from collections import Counter


"""Completes y = M*x + c conversions for channel numbers."""
def loadChannelEnergie(file):
    from spectrometryConfigurations import c, M

    return (M * channel) + c


"""Completes y = M*x + c conversions for channel numbers."""
def linear_channel_to_energy(channel):
    from spectrometryConfigurations import c, M

    return (M * channel) + c

"""Completes y = ax^2 + bx + c conversions for channel numbers."""
def parabolic_channel_to_energy(channel):
    from spectrometryConfigurations import a, b, c

    return (a * channel ** 2) + (b * channel) + c

"""Find most common occurrence in a list."""
def most_common(lst):
    return max(set(lst), key=lst.count)

"""Finds closest node to specified value in an array and returns the index."""
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

"""Finds closest X nodes to specified value in an array and returns the indices."""
def find_nearestX(array, value, x):
    idxs = []
    for i in range(x):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        idxs.append(idx)
        array = np.delete(array, idx)
    return idxs

"""Function which defines a gaussian curve."""
def gaussian(x, mu, sig, a):
    return a * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

"""Completes y = M*x + c conversions for channel numbers."""
def channel_to_energy(channel):
    return (read_calibration()[1] * channel) + read_calibration()[0]

"""Read config file to get multiplier and constant."""
def read_calibration():
    from spectrometryConfigurations import c, M
    return [c, M]

"""Corrects input spectrum for detecting efficiency at different energies based on scintillator material."""
def efficiency_correction(pathEnergy, array):
    efficiency_ref = pd.read_csv(pathEnergy+'Energy-Efficiency Relations_split.csv',delimiter=';')
    df = pd.DataFrame(efficiency_ref)
    corrected = np.zeros(len(array))

    from spectrometryConfigurations import scintillator, size

    for i in range(len(df["Scintillator"])):
        if scintillator == df["Scintillator"][i] and (np.abs(df["Thickness (mm)"][i] - size)) == 0:
            A1, t1, y0 = df["A1"][i], df["t1"][i], df["y0"][i]

    efficiency = A1 * (np.e ** (-(channel_to_energy(array) * 0.001) / t1)) + y0

    #corrected[i] = array[0] / efficiency
    corrected = array / efficiency
    now = str(datetime.today())
    todayDate=now.split(' ')[0]

    plt.plot(array,label="Raw Data")
    plt.plot(corrected,label="Corrected")
    plt.legend()
    pathOutput = pathEnergy.replace('input', 'output')
    pathFigure = pathOutput + todayDate
    if not os.path.exists(pathFigure):
        os.makedirs(pathFigure)
    fileFigure= pathFigure+ '/correctedVSNotCorrected.png'
    plt.savefig(fileFigure)
    plt.close()

    return corrected

"""Finds peaks based on Mariscotti's method (M.A. Mariscotti, Nucl. Instrum. Method 50 (1967) 309.)."""
def peak_finding(N, z, m, FWHM):
    # take second difference and standard deviation of spectrum
    S = second_diff(N, m, z)
    F = standard_dev_S(N, m, z)

    # define empty array same size as data
    signals = np.zeros(len(N))

    # get confidence factor from config
    from spectrometryConfigurations import confidence, intensity_threshold

    # send spectra from array to a list
    x = N.tolist()

    # smooth
    for i in range(z):
        x = average(x, m)

    # if the second difference is negative (indicates gaussian-like feature) and magnitude is greater than
    # the standard deviation multiplied by the confidence factor then can say a peak has been found,
    # if peak is greater than a % of spectrum's intensity, set index of peak centre to 1 in the signals array
    for i in range(len(N)):
        if abs(S[i]) > F[i] * confidence and S[i] < 0:
            try:
                if x[i] == max(x[i - FWHM:i + FWHM]) and x[i] >= (intensity_threshold / 100) * max(x):
                    signals[i] = 1
            except ValueError:
                pass

    return signals

"""Finds second difference of spectrum, discrete analogue to second derivative."""
def second_diff(N, m, z):
    # empty list with same size as data
    S = np.zeros(len(N))

    # loop over data points and for each take the second difference
    for i in range(len(N)):
        try:
            S[i] = N[i + 1] - 2 * N[i] + N[i - 1]
        # if at the ends of the array need to do different calculation to avoid error
        except IndexError:
            if i == 0:
                S[i] = N[i + 1] - 2 * N[i]
            else:
                S[i] = N[i - 1] - 2 * N[i]

    # smooth the second difference using Mariscotti's values
    for i in range(z):
        S = average(S, m)

    return S

"""Smooths input function by averaging neighboring values."""
def average(A, m):

    # define empty array same length as data
    B = np.zeros(len(A))

    # loop over data points and sum all points m either side of current one
    for i in range(len(A)):
        for j in range(i - m, i + m):
            try:
                B[i] += A[j]
            # if tries to go outside of boundaries of data, do nothing
            except IndexError:
                pass

    return B

"""Finds standard deviation of second difference."""
def standard_dev_S(N, m, z):
    # define empty array same size as data
    F = np.zeros(len(N))

    # loop over data points and take variance of the second difference
    for i in range(len(N)):
        try:
            F[i] = N[i + 1] + 4 * N[i] + N[i - 1]
        # need different calculation at ends of data array to avoid error
        except IndexError:
            if i == 0:
                F[i] = N[i + 1] + 4 * N[i]
            else:
                F[i] = N[i - 1] + 4 * N[i]

    # smoothing
    for i in range(z):
        F = average(F, m)

    # return square root of variance, std deviation
    return np.sqrt(F)

"""Removes straight line background from underneath the peaks."""
def remove_bg(array, locs, FWHM):

    correction = np.zeros(array.shape)
    for i in range(len(locs)):
        start, end = choose_region_around_index(array, locs, i, int(FWHM * 2))
        if start > 5 and end < len(array) - 6:
            # find mean count near start and end of slice
            before = np.mean(array[start - 5:start + 5])
            after = np.mean(array[end - 5:end + 5])
        elif start <= 5:
            before = np.mean(array[0:start + 5])
            after = np.mean(array[end - 5:end + 5])
        elif end >= len(array) - 6:
            before = np.mean(array[start - 5:start + 5])
            after = np.mean(array[end - 5:end])

        # define correction as straight line between mean count before and after peak
        correction_vals = np.linspace(before, after, end-start)
        correction[start:end] = correction_vals.astype(int)

        corrected = array - correction  # correct the array

    # make sure none of the array is negative
    if any(a <= 0 for a in corrected):
        for i in range(len(corrected)):
            if corrected[i] <= 0:
                corrected[i] = 0

    return corrected, correction

"""Calculate start and end indices for a slice around a given data point and given the size of the slice.
    Also accounts for index errors possible near the first and last index of the array."""
def choose_region_around_index(array, idxs, idx, val):
    try:
        # try to define start and end indices for slice
        start = idxs[idx] - val
        end = idxs[idx] + val
        # if one of the indices is out of bounds need to correct
        if start < 0 or end >= len(array):
            raise ValueError
    except ValueError:
        if idxs[idx] - val < 0:
            start = 0
            end = idxs[idx] + val    # if start is below 0 start slice at 0 to avoid index error
        elif idxs[idx] + val >= len(array):
            start = idxs[idx] - val
            end = len(array) - 1    # similar for end of array, end at the final index instead to avoid error

    return start, end

"""Reads lookup table and attempts to match peaks to isotopes."""
def source_lookup(pathFile, Es, r2s, peak_areas, corrected):
    try:
        from spectrometryConfigurations import time
    except ImportError:
        print("Time value not provided, CPS reported will be total counts")
        time = 1

    from spectrometryConfigurations import category, energy_window_size
    energy_ref = pd.read_csv(pathFile,delimiter=';')  # read energy reference
    df = pd.DataFrame(energy_ref)   # turn energy ref into pandas data frame
    d = {}

    # convert the peak areas to relative intensities and their errors as well as counts per second for each peak
    intensities = findintensity(peak_areas)
    Is = intensities[::2]
    I_errs = intensities[1::2]
    cps = peak_areas[::2] / time

    # generate empty lists of values
    sources = [0] * len(Es)
    isotopes = [0] * len(Es)
    absoluteIs = [0] * len(Is)
    ACFs = [0] * len(Is)

    # add all isotopes with a peak within 5% of calculated peak energy to a dataframe with unique name
    # for each identified peak this is then added to a dictionary for use when trying to identify which isotope
    # is responsible for each peak
    if category == "All":
        for i in range(len(Es)):
            i_list = []
            for j in range(len(df["Energy (keV)"])):
                if Es[i] - energy_window_size <= df["Energy (keV)"][j] <= Es[i] + energy_window_size:
                    i_list.append(df.loc[j])
            d["df" + str(i)] = pd.DataFrame.from_records(i_list, columns=df.columns)
    else:
        for i in range(len(Es)):
            i_list = []
            for j in range(len(df["Energy (keV)"])):
                if Es[i] - energy_window_size <= df["Energy (keV)"][j] <= Es[i] + energy_window_size \
                        and df["Category"][j] == category:
                    i_list.append(df.loc[j])
            d["df" + str(i)] = pd.DataFrame.from_records(i_list, columns=df.columns)

    # add all dataframes in the dicitonary together (one big list of all possible peaks for the identified peaks)
    dfx = pd.DataFrame(columns=df.columns)
    for i in range(len(Es)):
        dfx = pd.concat([dfx, d["df" + str(i)]], ignore_index=True)

    # try to identify the isotopes responsible for each peak
    for i in range(len(Es)):
        sources[i], isotopes[i], absoluteIs[i], ACFs[i] = find_optimum_source(d, i, Es, Is)

    # if multiple peaks and the spectrum is efficiency corrected, check the intensity ratios of the peaks
    if len(Es) > 1 and corrected:
        Es, Is, I_errs, absoluteIs, isotopes, sources, cps, r2s, ACFs, d, hidden_peaks = \
            check_ratios(Es, Is, absoluteIs, I_errs, isotopes, sources, cps, r2s, ACFs, d)

        # if hidden peaks found need to reassign isotopes to the peaks
        if hidden_peaks == True:
            for i in range(len(Es)):
                sources[i], isotopes[i], absoluteIs[i], ACFs[i] = find_optimum_source(d, i, Es, Is)

    # return the data for the output table
    identified_peaks = []
    for i in range(len(Es)):
        try:
            identified_peaks.append([sources[i], round(Es[i], 2),
                                    str(np.round(Is[i], 2)) + u"\u00B1" +
                                    str(np.round(I_errs[i], 2)),
                                    round(r2s[i], 4), round(cps[i], 2), ACFs[i]])
        except:
            print(I_errs[i])

    return identified_peaks

"""check peak intensity ratios and attempt to find hidden peaks"""
def check_ratios(Es, Is, absoluteIs, I_errs, isotopes, sources, cps, r2s, ACFs, d):
    # assume no hidden peaks initially
    hidden_peaks = False

    # for each peak calculate the ratio of intensities and compare to the expected ratio of intensities
    for i in range(len(Es)):
        for j in range(len(Es)):
            if isotopes[i] == isotopes[j] and i != j and isotopes[i] != 0:
                A = absoluteIs[i] / absoluteIs[j]
                B = Is[i] / Is[j]#

                # if one is found need to add the extra hidden peak to the lists of quantities
                if A - (A * 0.4) <= B <= A + (A * 0.4):
                    pass
                else:
                    hidden_peaks = True
                    if A > 1 and B > 1 or A < 1 and B < 1:
                        A_to_R = A / B
                        if A_to_R > 1:
                            Es = np.append(Es, Es[j])
                            Is = np.append(Is, Is[j] - (Is[j] / A_to_R))
                            Is[j] = Is[j] / A_to_R
                            d["df" + str(len(Es) - 1)] = d["df" + str(j)]
                            isotopes.append(0)
                            sources.append(0)
                            absoluteIs.append(0)
                            r2s = np.append(r2s, r2s[j])
                            cps = np.append(cps, cps[j])
                            I_errs = np.append(I_errs, I_errs[j])
                            ACFs.append(1)
                            break
                        else:
                            Es = np.append(Es, Es[i])
                            Is = np.append(Is, Is[i] - (Is[i] * A_to_R))
                            Is[i] = Is[i] * A_to_R
                            d["df" + str(len(Es) - 1)] = d["df" + str(i)]
                            isotopes.append(0)
                            sources.append(0)
                            absoluteIs.append(0)
                            r2s = np.append(r2s, r2s[i])
                            cps = np.append(cps, cps[i])
                            I_errs = np.append(I_errs, I_errs[i])
                            ACFs.append(1)
                            break
                    else:
                        pass
            break

    return Es, Is, absoluteIs, I_errs, isotopes, sources, cps, r2s, ACFs, d, hidden_peaks


"""Finds peak intensity relative to other peaks in spectra."""
def findintensity(areas):
    try:
        return 100 * (areas / (max(areas)))
    except ValueError:
        return []

"""Fits the peaks with a gaussian and returns the fit plus r^2 test statistic."""
def peak_fitting(array, locs, FWHM_est, removed):
    # define empty arrays for r^2 vals and gaussian fits
    new_locs, std_devs, r2s,  FWHMs, net_areas, gross_areas = [], [], [], [], [], []
    fits = np.zeros(array.shape)

    # cycle through peaks to fit each one
    for i in range(len(locs)):
        start, end = choose_region_around_index(array, locs, i, int(2 * FWHM_est))     # get start and end indices of slice

        data = array[start:end]     # take slice of spectrum

        if len(data) < 3:
            continue

        # x values to go into gaussian equation
        x = np.linspace(start, end, end-start)

        # estimated parameters
        mean = locs[i]
        sigma = FWHM_est / 2.355
        amp = max(data)

        # calculate gaussian fit parameters from given estimated parameters
        try:
            popt, pcov = curve_fit(gaussian, x, data, p0=[mean, sigma, amp])
            new_loc = popt[0]
            std_dev = abs(popt[1])
            FWHM = abs(popt[1]) * 2.355

            #new_start, new_end = choose_region_around_index(array, [np.around(new_loc).astype(int)], 0, int(2 * FWHM))
            #new_data = array[new_start:new_end]

            # create gaussian fit
            x_fit = np.linspace(start, end, end - start)
            y_fit = gaussian(x_fit, *popt)
            fits[start:end] = y_fit

            r2 = coefficient_of_determination(data, y_fit, len(data))
        except:
            r2 = 0

        fits[start:end] += removed[start:end]  # background adjustment to fits to match spectrum

        # if the fit is bad, reject the peak
        if not r2 < 0.1:
            new_locs.append(new_loc)
            r2s.append(r2)
            FWHMs.append(FWHM)
            net_areas.append(np.sum(y_fit))
            net_areas.append(np.sqrt(np.sum(y_fit)))
            gross_areas.append(np.sum(y_fit + removed[start:end]))
            gross_areas.append(np.sqrt(np.sum(y_fit + removed[start:end])))
        else:
            fits[start:end] = 0
            try:
                fits[start:end] = 0
            except:
                pass

    return np.array(new_locs), np.array(r2s), fits, np.array(net_areas), np.array(gross_areas), np.array(FWHMs)

"""Calculate the coefficient of determination (r^2 value) for a given data set and fit."""
def coefficient_of_determination(y, y_fit, n):
    # calculate residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)

    # calculate total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2) + 1E-9

    # calculate r-squared
    r2 = 1 - (ss_res / ss_tot)

    # number of variables in the data
    k = 2

    # return adjusted r^2
    return 1 - (((1 - r2) * (n - 1)) / (n - k - 1))

"""Finds the most likely source causing a peak based on different criteria."""
def find_optimum_source(d, i, Es, Is):
    from spectrometryConfigurations import common_isotopes

    # if no sources in dataframe for this peak then skip it
    dfi = d["df" + str(i)]
    if len(dfi["Isotope"]) == 0:
        return 0, 0, 0, 0

    # empty lists to be added to
    sources, isotope_list, absoluteIs = [], [], []

    # first check if multiple peaks could be same source (are there same isotopes in multiple of the distinct dfs)
    if len(Es) > 1:
        source_list = []
        for j in range(len(Es)):
            if i != j:
                dfj = d["df" + str(j)]
                list_i = set(dfi["Isotope"])
                list_j = set(dfj["Isotope"])
                for l in list_i:
                    for k in list_j:
                        if l == k:
                            source_list.append(l)

        # try assigning the most common occurrence in list of isotopes occurring in multiple of the distinct dfs
        try:
            isotope_list.append(most_common(source_list))
            for j in range(len(dfi["Isotope"])):
                if dfi["Isotope"][j] == isotope_list[0]:
                    sources.append(dfi["Isotope"][j] + " (" + dfi["Sub-category"][j] + ")")
                    absoluteIs.append(dfi["Absolute emission probability (%)"][j])
        except ValueError:
            pass

    # then see if any of the common isotopes could be attributed
    if any(a in common_isotopes for a in dfi["Isotope"]):
        temp_Es, temp_Is, temp_isotopes, temp_sources = [], [], [], []
        for j in range(len(dfi["Energy (keV)"])):
            if dfi["Isotope"][j] in common_isotopes:
                temp_Es.append(dfi["Energy (keV)"][j])
                temp_Is.append(dfi["Absolute emission probability (%)"][j])
                temp_isotopes.append(dfi["Isotope"][j])
                temp_sources.append(dfi["Isotope"][j] + " (" + d["df" + str(i)]["Sub-category"][j] + ")")

        # try assign common isotope with peak closest in energy to the identified peak
        if len(temp_Es) > 0:
            index = find_nearest(temp_Es, Es[i])
            sources.append(temp_sources[index])
            isotope_list.append(temp_isotopes[index])
            absoluteIs.append(temp_Is[index])

    # then try to find isotope with peak closest to the found peak's energy
    if len(dfi["Energy (keV)"]) > 3:
        closest_Es = find_nearestX(dfi["Energy (keV)"], Es[i], 3)
        closest_Is = find_nearestX(dfi["Absolute emission probability (%)"], Is[i], 3)
        idx = None
        for a in closest_Es:
            if a in closest_Is:
                idx = a

        if idx == None:
            idx = closest_Es[0]

        sources.append(dfi["Isotope"][idx] + " (" + dfi["Sub-category"][idx] + ")")
        isotope_list.append(dfi["Isotope"][idx])
        absoluteIs.append(dfi["Absolute emission probability (%)"][idx])
    else:
        idx = find_nearest(dfi["Energy (keV)"], Es[i])
        sources.append(dfi["Isotope"][idx] + " (" + dfi["Sub-category"][idx] + ")")
        isotope_list.append(dfi["Isotope"][idx])
        absoluteIs.append(dfi["Absolute emission probability (%)"][idx])

    # how many of the above methods agree on the isotope determines the confidence factor
    arbitrary_confidence_factor = 1
    if len(isotope_list) == 3:
        if len(set(isotope_list)) == 1:
            arbitrary_confidence_factor = 3
        elif len(set(isotope_list)) == 2:
            arbitrary_confidence_factor = 2
    elif len(isotope_list) == 2:
        if len(set(isotope_list)) == 1:
            arbitrary_confidence_factor = 2
    else:
        arbitrary_confidence_factor = 1

    if len(Counter(isotope_list))>1:
        # return the first one in the list as this is probably the best answer
        return sources[np.argmax(absoluteIs)], isotope_list[isotope_list.index(sources[np.argmax(absoluteIs)].split(' ')[0])],\
               absoluteIs[np.argmax(absoluteIs)], arbitrary_confidence_factor
    else:
        # return the first one in the list as this is probably the best answer
        return sources[0], isotope_list[0], absoluteIs[0], arbitrary_confidence_factor