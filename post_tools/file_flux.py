import numpy as np
index=0
list_total=[]
list_electron=[]
list_positron=[]
with open("outp","r") as f:
    line=f.readline()
    while line:
        if line[:9]==' cell (10':
            f.readline()
            index=index+1                                                                                                                    
            for i in range(201):
                energy,f4,error=f.readline().split()
                if((index%3)==1):
                    list_total.append(f4)
                elif((index%3)==2):
                    list_electron.append(f4)
                else:
                    list_positron.append(f4)
        line=f.readline()

print(list_total)

