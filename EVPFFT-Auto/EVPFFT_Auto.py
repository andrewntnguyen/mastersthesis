# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 17:22:44 2021

@author: thean
"""
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
from scipy import interpolate as inter

#---------------Parameters Pickup-------------------------#
fpara = "parameters.txt" #Parameter file
with open(fpara, 'r') as fpar:
    para = fpar.readlines()
load_dir = para[1][2]

# Sample Voxel Size 32 64 128#
para3 = para[3].split("\t")
para5 = para[5].split("\t")
para7 = para[7].split("\t")
SampleSize = para7[1]
nsteps = para7[0]
# Sample ID Name#
SampleID = para7[2]

fname = "cupl2_"+SampleSize+"cubed.sx" # CUPL2 Initial for each sample size has the calibrated Tau0, TauS, Theta0, ThetaS
with open(fname, 'r') as file:
    # read a list of lines into data
    data = file.readlines()

finp = "fft_"+load_dir+"_init.in" #fft parameter file
with open(finp, 'r') as ffft:
    fft = ffft.readlines()

fdim = "fft_init.dim" #dimension file
with open(fdim, 'r') as filedim:
    dim = filedim.readlines()

mode = int(para[1][0]) # Mode = 0: Range and # of iterations; Mode = 1: Tabular Data; Mode = 2: Multifft

#---------Rewrite Input files for EVPFFT-----------------------------------------------------#
if mode == 0 or mode == 1:
    fft[1] = str(int(para7[1])**3)+'\t\t\tnumber of Fourier points (should be npts1*npts2*npts3)\n' #Number of Fourier Points
    fft[4] = str(para7[3]+'\n') #FFT File name
    fft[31] = str(nsteps)+'\t\tnsteps\n' #N steps 
    fft[39] = str('1 ')+str(para7[0])+'            IWFIELDS,IWSTEP\n' #IWStep
    dim[4]= "      PARAMETER(NPTS1="+str(SampleSize)+",NPTS2="+str(SampleSize)+",NPTS3="+str(SampleSize)+")\n" 
    fftout = "fft_"+load_dir+".in"
    fftdim = "fft.dim"
    
    with open(fftout,'w') as file:
        file.writelines( fft )
    with open(fftdim, 'w') as file:
        file.writelines( dim )
        
#Mode0 --- Varried Number of Runs
if mode == 0:
    val = para[12].split("\t")
    VAR=0
    for i in range(int(para5[0])):
        Tau0 = float(val[0])+i*float(para3[0]) #Tau0 varies
        for j in range(int(para5[1])):
            TauS = float(val[1])+j*float(para3[1]) #TauS varies
            for k in range(int(para5[2])):
                Theta0 = float(val[2])+k*float(para3[2]) #Theta0 varies
                for l in range(int(para5[3])):
                    Theta1 = float(val[3])+l*float(para3[3]) #Theta1 varies
                    data[8] = '  '+str(Tau0)+' '+str(Tau0)+' '+str(TauS-Tau0)+' '+str(Theta0)+' '+str(Theta1)+'\ttau0xf,tau0xb,tau1x,thet0,thet1\n'
                    fout1 = "fft_"+load_dir+str(i+1)+".in"
                    fout2 = "cupl2_"+str(i+1)+".sx"
                    with open(fout1, 'w') as file:
                        file.writelines( fft )
                    with open(fout2, 'w') as file:
                        file.writelines( data )
                    VAR+=1
                    
#Mode1 --- Tabular Data for Tau0, TauS, Theta0, ThetaS
if mode == 1:
    d = np.loadtxt(fpara, dtype=str, delimiter='\t', skiprows=12)
    VAR = len(d)-1
    for i in range(VAR):
        if d[i][0] != "Stop":
            Tau0 = float(d[i][0])
            TauS = float(d[i][1])
            Theta0 = float(d[i][2])
            Theta1 = float(d[i][3])
            data[8] = '  '+str(Tau0)+' '+str(Tau0)+' '+str(TauS-Tau0)+' '+str(Theta0)+' '+str(Theta1)+'\ttau0xf,tau0xb,tau1x,thet0,thet1\n'
            fout1 = "fft_"+load_dir+str(i+1)+".in"
            fout2 = "cupl2_"+str(i+1)+".sx"
            with open(fout1, 'w') as file:
                file.writelines( fft )
            with open(fout2, 'w') as file:
                file.writelines( data )

#Mode2 --- Multiple FFTs
if mode == 2:
    fftlist = "fftlist.txt" #List of fft file names
    d = np.loadtxt(fftlist, dtype=str, delimiter='\t', skiprows=6)
    VAR = len(d)-1
    for i in range(VAR):
        fft[1] = str(int(para7[1])**3)+'\t\t\tnumber of Fourier points (should be npts1*npts2*npts3)\n' #Number of Fourier Points
        fft[4] = str(d[0]+d[i+1]+'-fft.txt\n') #Name of FFT File
        fft[31] = str(para7[0])+'\t\tnsteps\n' #N steps 
        fft[39] = str('1 ')+str(para7[0])+'            IWFIELDS,IWSTEP\n' #IWStep
        dim[4]= "      PARAMETER(NPTS1="+str(SampleSize)+",NPTS2="+str(SampleSize)+",NPTS3="+str(SampleSize)+")\n" 

        fftdim = "fft.dim"
        fout1 = "fft_"+load_dir+str(i+1)+".in"
        fout2 = "cupl2_"+str(i+1)+".sx"
        with open(fout1, 'w') as file:
            file.writelines( fft )
        with open(fftdim, 'w') as file:
            file.writelines( dim )
        with open(fout2, 'w') as file:
            file.writelines( data )
                
#--------------------Run Simulation-----------------------------------------------#
Directory1 = str(SampleSize)+str(SampleID)+"_"+load_dir
os.system('mkdir '+Directory1)
for i in range(VAR):
    if mode == 2:
        #Deg = [d[i+1][j:j+2] for j in range(0,len(d[i+1]),2)]
        Label = d[i+1]
        SubDirect1 = Directory1+"/"+Label
        #SubDirect1 = Directory1+"/"+Deg[0]+Deg[1]+Deg[2]
        SubDirect1A = Directory1+"-"+Label
        #SubDirect1A = Directory1+"-"+Deg[0]+Deg[1]+Deg[2]
    if mode == 1 or mode == 0:
        SubDirect1 = Directory1+"/"+d[i][4]
        SubDirect1A = Directory1+"-"+d[i][4]
    os.system('mkdir '+SubDirect1)
    
    loadtag =load_dir+str(i+1)         #load direction + Iteration number
    # fft_X/Y/Z.in
    os.system('cp fft_'+loadtag+".in "+SubDirect1)
    os.system('mv fft_'+loadtag+".in fft_"+load_dir+".in")
    
    # fft.dim
    os.system('cp fft.dim '+SubDirect1)
    
    # cupl2.sx
    os.system('cp cupl2_'+str(i+1)+'.sx '+SubDirect1)
    os.system('mv cupl2_'+str(i+1)+'.sx cupl2.sx')
    
    # cuel2.sx
    os.system('cp cuel2.sx '+SubDirect1)

    # EVPFFT
    print(' Start of EVPFFT for '+SubDirect1+'-------------------------------------------------\n')
    
    os.system('gfortran evpnew10_'+load_dir+'2.for')
    os.system('./a.out')
    
    # str_str file save
    os.system('mv str_str.out str_str_'+Directory1+'_'+str(i+1)+'.out')
    os.system('cp str_str_'+Directory1+'_'+str(i+1)+'.out '+SubDirect1)
    
    # Field files save
    os.system('cp dfield.out '+SubDirect1)
    os.system('cp efield.out '+SubDirect1)
    os.system('cp elfield.out '+SubDirect1)
    os.system('cp sfield.out '+SubDirect1)
    
    # cupl2.sx for Plots
    os.system('mv cupl2.sx cupl2_'+str(i+1)+'.sx')
    print('\nEnd of EVPFFT for '+SubDirect1+'---------------------------------------------------\n')


#------------------Post process plots and RMS-------------------------------------------#
#----------Experimental Data----------#
Eexp = list()
Sexp = list()
theo = para[9].split("\t")
theo_name = theo[0]         # Name of experimental data csv file in Parameters file
with open(theo_name,encoding='utf-8-sig') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    for row in readCSV:
        Eexp.append(float(row[0]))
        Sexp.append(float(row[1]))
f = inter.interp1d(Eexp,Sexp)
plt.plot(Eexp,Sexp,label='Experimental')

#---------RambergOsgood---------#
#RamOsPara = para[9].split("\t")
#Young = float(RamOsPara[0])
#Sty = float(RamOsPara[1])
#N = float(RamOsPara[2])
#Ult_S = float(RamOsPara[3])
#RamOsS = []
#RamOsE = []

#for i in range(0,int(nsteps)):
#    Stress_tmp = Ult_S/int(nsteps)*i
#    Strain_tmp = Stress_tmp / Young + 0.002*(Stress_tmp / Sty)**(1/N)
#    RamOsE.append(Strain_tmp)
#    RamOsS.append(Stress_tmp)
#plt.plot(RamOsE,RamOsS,label='Ramberg-Osgood')

#-------------Plotting + PercentError&RMS-----------------#
#Directory1 = str(SampleSize)+str(SampleID)+"_"+load_dir
#os.system('mv str_str.out str_str_'+Directory1+'_'+str(i+1)+'.out')
file = open("RMS_"+SubDirect1A+'.txt',"w")
if load_dir == 'Z': #S33
    Start = 228
    End = 241
if load_dir == 'Y': #S22
    Start = 216
    End = 229
if load_dir == 'X': #S11
    Start = 204
    End = 217
for m in range(VAR):
    fname = 'str_str_'+Directory1+'_'+str(m+1)+'.out'
    d = np.loadtxt(fname, dtype=str, delimiter='\t', skiprows=1)
    stress=''
    strain=''
    B=0
    val_S = [0]
    val_E = [0]
    PerErr = []
    for i in range(0,int(nsteps)): # Go through each rows 

        for j in range(Start,End): #Range of Uniaxial Stress characters
            if j % 12 !=0:
                if d[i][j] !=' ':
                    stress = stress+d[i][j]
                    #print(stress)
                else:
                    B+=1
                    val_S.append(float(stress))
                    stress=''
        for k in range(36,49): #Range of E33's characters
            if k % 12 !=0:
                if d[i][k] !=' ':
                    strain = strain+d[i][k]
                    #print(strain)
                else:
                    val_E.append(float(strain))
                    strain = ''
        #print(val_S)
        #print(val_E)
    with open('cupl2_'+str(m+1)+'.sx','r') as fcupl2:
        cupl = fcupl2.readlines()
    ParaTauTheta = cupl[8]
    
    if mode == 0:
        A = ParaTauTheta.split("\t")
        B = A[0].split(" ")
        plt.plot(val_E,val_S,label=B[2]+' '+B[3]+' '+B[4]+' '+B[5])
    if mode == 1:    
        d = np.loadtxt(fpara, dtype=str, delimiter='\t', skiprows=12)
        plt.plot(val_E,val_S,label=d[m][4])
    if mode == 2:
        deg = np.loadtxt(fftlist, dtype=str, delimiter='\t',skiprows=7)
        plt.plot(val_E,val_S,label=deg[m])  
    
    SumSquare = 0
    PerErr=[]
    for i in range(1,int(nsteps)+1): 
        #val_expected = val_S[i] / Young + 0.002*(val_S[i] / Sty)**(1/N)
        val_expected = f(val_E[i])
        PerErr.append(abs((val_S[i]-val_expected)/val_expected)*100)
        SumSquare += (abs((val_S[i]-val_expected)/val_expected)*100)**2
    RMS = np.sqrt(SumSquare/int(nsteps))   
    print('Average percent error of set #'+str(m)+': '+ str(sum(PerErr)/len(PerErr))+' %')
    print('Root Mean Square of try #'+str(m)+': '+ str(RMS))

    file.write(ParaTauTheta)
    file.write('Avg % error of Iteration#'+str(m)+': '+ str(sum(PerErr)/len(PerErr))+'%\nRMS of Iteration#'+str(m)+': '+ str(RMS)+'\n\n') 
    os.system('rm cupl2_'+str(m+1)+'.sx')
    os.system('rm '+fname)
file.close()



#-----------ShowPlot------------#
plt.xlabel('Strain')
plt.ylabel('Stress')
plt.title(Directory1+' Parametric Study')
plt.legend()
plt.savefig('str_str_'+SubDirect1A+'.png')
plt.show()

SubDirect2 = Directory1+'/'+"ArchivePlot"
SubDirect3 = Directory1+'/'+"ArchiveRMS"
os.system('mkdir '+SubDirect2)
os.system('mkdir '+SubDirect3)
os.system('cp str_str_'+SubDirect1A+'.png '+SubDirect2) #Save stress strain plot to 
os.system('rm str_str_'+SubDirect1A+'.png')
os.system('cp RMS_'+SubDirect1A+'.txt '+SubDirect3)
os.system('rm RMS_'+SubDirect1A+'.txt')
os.system('rm fft_'+load_dir+'.in fft.dim')
os.system('rm dfield.out efield.out elfield.out sfield.out a.out conv.out vm.out tex.out err.out')