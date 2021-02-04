import matplotlib.pyplot as plt
import csv

EZ = list()
SAZ = list()
SBZ = list()
SCZ = list()
SDZ = list()
SEZ = list()
EExZ = list()
SExZ = list()
with open('515-2_75x_CTFCuts.csv',encoding='utf-8-sig') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    for row in readCSV:
        print(row)
        EZ.append(float(row[0]))
        SAZ.append(float(row[1]))
        SBZ.append(float(row[2]))
        EExZ.append(float(row[3]))
        SExZ.append(float(row[4]))
        #ECZ.append(float(row[4]))
        #SCZ.append(float(row[5]))
        #EDZ.append(float(row[6]))
        #SDZ.append(float(row[7]))
        #EEZ.append(float(row[8]))
        #SEZ.append(float(row[9]))

plt.plot(EZ,SAZ,label='Sample A')
plt.plot(EZ,SBZ,label='Sample B')
#plt.plot(ECZ,SCZ,label='Load Z')
plt.plot(EExZ,SExZ,label='Ji 2020 experimental')
#plt.plot(EEZ,SEZ,label='32E')

plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title('CTF Cuts')
plt.legend()
plt.show()
