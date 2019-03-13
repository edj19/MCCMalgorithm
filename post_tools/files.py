import os
import numpy as np

def eachFile(filepath):
	j=0
	pathDir = os.listdir(filepath)   #获取当前路径下的文件名，返回list
	for s in pathDir:
		newDir=os.path.join(filepath,s)   #将文件名写入到当前文件路径后面
		if os.path.isfile(newDir): #如果是文件
		    if os.path.splitext(newDir)[1]==".txt":  #判断是否是txt
		        readFile(newDir)
		        j=j+1
		        print(j)
		        pass
		    else:
		    	break
			# eachFile(newDir)



def readFile(filepath):
	# f_energy = open("energy.txt","w")
	# i=0
	index=0                  #控制数据存入不同的list
	# list_total=[]
	# list_electron=[]
	# list_positron=[]
	with open(filepath,"r") as f:
		line=f.readline()
# 		print(line)
		while line:
			if line[:9]==' cell  10':     #根据关键词抽取数据
				f.readline()
				index=index+1
				for i in range(126):       #抽取的数据格式
					energy,f4,error=f.readline().split()
					if(index==1):
						list_total.append(f4)
					elif(index==2):
						list_electron.append(f4)
					else:
						list_positron.append(f4)
					# f_energy.write(str(f4)+' ')
					# f_energy.write('\n')
				# i=i+1
				# energy,f4,error=f.readline().split()
				# # flux4.append(f4)
				# f_energy.write(str(f4)+' ')
				# f_energy.write('\n')
				# if (i%125==0):
				# 	break
				# else:
				# 	fout.write(str(f4)+' ')
				# # f.readline()
			line=f.readline()
# 		print(list_electron)
	# f_energy.close()
	# print(i)




def main():
	global list_total,list_electron,list_positron      #定义全局变量，可以将所有数据都存入list中
	fp=r'F:\\MCwork\\MCCM\\scripts\\filesworks'
	os.chdir(fp)
	eachFile(fp)
	output =open("flux.txt",'w')    #将list存入相应的文件中，便于后期处理
	print(len(list_electron))
	listdata_total=list(np.reshape(list_total,(15,126)).T)    #改变数组维度，存储
	for i in range(126):                      #数据读入相应文件的第一种方法，第一篇博客有介绍
		for j in range(15):
			output.write(listdata_total[i][j]+' ')
			output.write('\t')
		output.write('\n')
	# ftest=open("test.txt","w")
	# list1 = [[1,2,3],[4,5,6]]
	# #numpy
	# np.savetxt('a.txt',list1)
	output.close()

	#pandas
	# df
	# ftest.write('The test one txt file\n')
	# ftest.write(str(list1))
	# ftest.close()




if __name__ == '__main__':
	list_electron=[]
	list_total=[]
	list_positron=[]
	main()