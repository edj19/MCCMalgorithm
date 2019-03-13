#！/usr/bin/python
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pylab as plt
# import pandas as pd

def dpa_mccm():
	dpa = []
	depth = []
	data_Ti = np.loadtxt("photonvec_Ti1000.txt")
	x = data_Ti[:,1]
	data_O = np.loadtxt("photonvec_O1000.txt")
	y = data_O[:,1]
	flux = np.loadtxt("energy_1000.txt")
	for i in range(15):
		number = flux[:,i]
		value = (sum(np.multiply(x,number))/3+sum(np.multiply(y,number)*2/3))*0.01
		# value = (sum(np.multiply(y,number))*2/3)*0.01
		# value = (sum(np.multiply(x,number))/3+sum(np.multiply(y,number)*2/3))*0.01
		dpa.append(value)
	for i in range(15):
		delta = (i+1)*0.1
		depth.append(delta)
	with open("result.txt",'w') as f:
		for i in range(len(dpa)):
			f.write(str(dpa[i]))
			f.write("\n")
	print(dpa)
	plt.plot(depth,dpa,'bo--')
	plt.show()


if __name__ == '__main__':
	dpa_mccm()
