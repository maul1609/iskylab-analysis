from netCDF4 import Dataset
import numpy as np
import numpy as np

def readData(fileName,modelStr="Model-Exp005"):
	nc=Dataset(fileName)

	data1=dict()
	data1 = {modelStr : \
		{"time" : nc['time'][:], "p" : nc['p'][:], \
		"t" : nc['t'][:], "rh": nc['rh'][:], \
		"ndrop" : nc['ndrop'][:], "ql" : nc['ql'][:], \
		"deff" : nc['deff'][:], \
		'rhoa': nc['p'][:]/nc['t'][:]/(8.314/29e-3)}}
		
	nc.close()
	return data1

if __name__== "__main__":
	data1=readData()

