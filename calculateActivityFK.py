import numpy as np



## zn-65

energy= 1115.52 #keV

pol1=1.941
pol2=-1.521*np.log(energy)
pol3=0.1139*np.log(energy)**2
pol4=-0.005429*np.log(energy)**3

attenuationCorrectedEcfficieny=np.exp(pol1+pol2+pol3+pol4)

halflifetime=21116000
#attenuationCorrectedEcfficieny = 1.7601*10**-3
volumn=1
netArea=1.682*10**6


time=54000#sec
part1=halflifetime /(np.log(2)*time)
part2= 1*(np.log(2)*time/halflifetime)
Kc = part1*part2
tw=950400


Kw= np.exp(-np.log(2)*tw/halflifetime)



#1 Mikrocurie [µCi] = 37 000 Becquerel [bq]


# currie
conversionFactor=37000

#bq
conversionFactor=1

emissionProbability=0.5075
c = netArea/(volumn*attenuationCorrectedEcfficieny*emissionProbability*time*Kw*Kc*conversionFactor)


print(c)