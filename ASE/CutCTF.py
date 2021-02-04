# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:06:35 2021

@author: thean
"""
ctfin = "515-2_75x_Full.ctf"

with open(ctfin,'r') as ctf:
    datain = ctf.readlines()

xmin = 384
xmax = 512
ymin = 256
ymax = 384
xcells = 512
ycells = 384
fout = open('515-2_75x_L.ctf','w+')
for i in range(0,15):
    if i == 4:
        fout.write('XCells\t'+str(xmax-xmin)+'\n')
    elif  i == 5:
        fout.write('YCells\t'+str(ymax-ymin)+'\n')
    else:
        fout.write(datain[i])
count = 14
for y in range(1,ycells+1):
    for x in range(1,xcells+1):
        count+=1
        if x > xmin and x <= xmax:
            if y > ymin and y <= ymax:
                #print(x,' ',y,' ',count)
                fout.write(datain[count])
                #print(x,y)
fout.close()
