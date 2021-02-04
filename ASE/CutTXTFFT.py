min_X = 256
max_X = 384
min_Y = 0
max_Y = 128
min_Z = 128
max_Z = 256
listb = []
listd = []
with open('Synthetic_515-2_75x_128F_G.txt','x') as f_out:
    with open('Synthetic_Rolled_F_fft.txt', 'r') as file:
        for row in file:
            Angle1, Angle2, Angle3, X, Y, Z, GrainID, H = row.split()
            if int(X) > min_X and int(X) <= max_X:
                if int(Y) > min_Y and int(Y) <= max_Y:
                    if int(Z) > min_Z and int(Z) <= max_Z:
                        count_X = int(X)-min_X
                        count_Y = int(Y)-min_Y
                        count_Z = int(Z)-min_Z
                        f_out.write(str(Angle1)+' '+str(Angle2)+' '+str(Angle3)+' '+str(count_X)+' '+str(count_Y)+' '+str(count_Z)+' '+str(GrainID)+' '+str(H)+'\n')
print('done')                    
                    
        
