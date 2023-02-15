import numpy as np




energy= 661.64 #keV

pol1=-126.4
pol2=77.16*np.log(energy)
pol3=-18.04*np.log(energy)**2
pol4=1.843*np.log(energy)**3
pol5=-0.07027*np.log(energy)**4


attenuationCorrectedEcfficieny=np.exp(pol1+pol2+pol3+pol4+pol5)

attenuationCorrectedEcfficieny = 1.7601*10**-3
volumn=1
netArea=9384.9
time=4020#sec
Kw=0.999998537
Kc=0.99540311


tw=6328800
halflifetime=9.521*10**8
Kw= np.exp(-np.log(2)*tw/halflifetime)

#1 Mikrocurie [µCi] = 37 000 Becquerel [bq]

part1=halflifetime /(np.log(2)*time)
part2= 1*(np.log(2)*time/halflifetime)
Kc = part1*part2
conversionFactor=37000
emissionProbability=0.8512
c = netArea/(volumn*attenuationCorrectedEcfficieny*emissionProbability*time*Kw*Kc*conversionFactor)
print(c)
x=0

