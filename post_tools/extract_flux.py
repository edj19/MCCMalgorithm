import os
import numpy as np
def extractflux():
    index=0
    list_total=[]
    list_electron=[]
    list_positron=[]
    with open("outp.txt", "r") as f:
        line=f.readline()
        while line:
            if line[:9]==' cell  10':
                f.readline()
                for i in range(126):
                    energy,f4,error=f.readline().split( )
                    list_total.append(f4)
            line=f.readline()
        # if line[:12]=='          Th':
        #     print(1)
        # else:
        #     print(0)
        print(line[:12])
        # while line:
        #     # if line[:9]==' cell  10':
        #     #     f.readline()
        #     #     index=index+1
        #     #     for i in range(126):
        #     #         energy,f4,error = f.readline().split()
        #     #         if(index==1):
        #     #             list_total.append(f4)
        #     #         elif(index==2):
        #     #             list_electron.append(f4)
        #     #         else:
        #     #             list_positron.append(f4)
        #     line=f.readline
        # print(list_total)

def main():
    # filename='outq'
    extractflux()
                            
if __name__ == '__main__':
    main()