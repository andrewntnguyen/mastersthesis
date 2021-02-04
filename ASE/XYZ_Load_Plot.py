import matplotlib.pyplot as plt
import pandas as pd
import csv
E = list()
S11 = list()
S22 = list()
S33 = list()
    
with open('515-2_75x_Synthetic-ZLoad.csv',encoding='utf-8-sig') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    for row in readCSV:
        print(row)
        E.append(float(row[0]))
        S11.append(float(row[1]))
        S22.append(float(row[2]))
        S33.append(float(row[3]))

plt.plot(E,S11,label='SSE1')
plt.plot(E,S22,label='SSE2')
plt.plot(E,S33,label='SSE3')

plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title('515-2_75x_Synthetic Z Load')
plt.legend()
plt.show()

