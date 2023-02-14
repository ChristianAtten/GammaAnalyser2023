#import mainScript
import analyzingFunctions
import matplotlib.pyplot as plt
import numpy as np

def most_common(lst):
    """Find most common occurrence in a list."""
    return max(set(lst), key=lst.count)


def plotOutcome(pathOutputData, file, y, df1, fits, FWHMs, net_areas, gross_areas):

    try:
        x = analyzingFunctions.linear_channel_to_energy(np.arange(len(y)))
        #x=0
        # new
        #testX = data_help.iloc[new_locs, 1]
        #testX2 = testX.T.to_numpy()
        #x = testX2
    except ImportError:
        x = analyzingFunctions.parabolic_channel_to_energy(np.arange(len(y)))

    # check discriminators are within the length of the data
    from config_sam import LLD, HLD
    if LLD < 0:
        LLD = 0
    if HLD > len(y):
        HLD = len(y)

    # make sure nothing is plot over the spectrum if there is no fit there
    for i in range(len(fits)):
        if fits[i] <= 0.1:
            fits[i] = float('nan')

    # print(f'time elapsed final: {time.time()-start}')

    # plot spectrum and the fits over it
    isotopes = df1["Isotope"]

    try:
        from config_sam import M
        df1["FWHM (keV)"] = np.around(FWHMs * M, 2)
    except ImportError:
        from config_sam import a, b
        df1["FWHM (keV)"] = np.around((a * FWHMs ** 2) + (FWHMs * b), 2)

    df1["Net Area"] = np.around(net_areas[::2])
    uncerts = []
    for i in range(len(gross_areas[1::2])):
        try:
            uncerts.append(np.around(gross_areas[1::2][i] / df1[r"$\bar{R}^2$ of fit"][i], 2))
        except KeyError:
            pass
    df1["Net Area Uncertainty"] = uncerts
    fileName=file.split('/')[-1]
    df1.to_csv(pathOutputData + fileName[:-4] + "_peaks.csv", index=False)

    df1 = df1.drop(columns=["ACF", "CPS"])

    if len(isotopes) <= 10:
        fig, axs = plt.subplots(2, 1)
        axs[0].plot(x[LLD:HLD], y[LLD:HLD])
        axs[0].plot(x[LLD:HLD], fits[LLD:HLD])
        # axs[0].set_title("Energy (keV)", fontsize=12)
        #axs[0].set_yscale('log')
        axs[0].set_ylabel("Counts (Total)", fontsize=14)
        axs[0].set_xlabel("Energy (keV)", fontsize=14)

        plt.figure()
        plt.plot(x[LLD:HLD], y[LLD:HLD], label="Raw Data")
        plt.plot(x[LLD:HLD], fits[LLD:HLD], label="Identified Peaks")
        plt.legend(loc="upper right")

        #plt.set_ylabel("Counts (Total)", fontsize=14)
        #plt.set_xlabel("Energy (keV)", fontsize=14)
        x = 0

        try:
            axs[1].axis('tight')
            axs[1].axis('off')
            table = axs[1].table(cellText=df1.values,
                                 rowLabels=None,
                                 colLabels=df1.columns,
                                 loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(14)
        except:
            pass
    else:
        plt.plot(x[LLD:HLD], y[LLD:HLD])
        plt.plot(x[LLD:HLD], fits[LLD:HLD])
        plt.yscale('log')
        plt.text(0.2 * max(x[LLD:HLD]), 0.9 * max(y[LLD:HLD]), str(len(isotopes)) +
                 " peaks found\nMost common isotope is " + str(most_common(isotopes.tolist())), fontsize=12)
        plt.ylabel("Counts (Total)", fontsize=12)
        plt.xlabel("Energy (keV)", fontsize=12)



    plt.show()
    fileFigure = pathOutputData + '/totalResult.png'
    plt.savefig(fileFigure)
