import matplotlib.pyplot as plt
import numpy as np
import readMeteoCPC
import readOPC_Merged
import readPNSD_Mrg
import svp
from scipy.interpolate import interp1d

cloud_times=dict()
cloud_times['Exp006']={'time' :np.array([14,20])*60.}
cloud_times['Exp007']={'time' : np.array([18,21])*60.}
cloud_times['Exp008']={'time' : np.array([15,18])*60.}
cloud_times['Exp028']={'time' : np.array([11,13])*60.}

def calc_average(cloud_times,meteo,opc):
	ind,=np.where((cloud_times['time'][0]<=meteo['Time']) &  \
		(cloud_times['time'][1]>meteo['Time']))
	cdnc=np.mean(opc['ndrop'][ind])
	
	return cdnc

if 'data1' in locals():
	pass
else:
	data1=dict()
	
"""
	Meteo OPC
"""
for i in range(len(readMeteoCPC.metStr)):
	data2=readMeteoCPC.readData(readThis=i,metStr=readMeteoCPC.metStr[i])
	data1[readMeteoCPC.metStr[i]]=data2[readMeteoCPC.metStr[i]].copy()
"""
	OPC data
"""
for i in range(len(readOPC_Merged.opcStr)):
	data2=readOPC_Merged.readData(readThis=i,opcStr=readOPC_Merged.opcStr[i])
	data1[readOPC_Merged.opcStr[i]]=data2[readOPC_Merged.opcStr[i]].copy()

"""
	NPSD data
"""
for i in range(len(readPNSD_Mrg.npsdStr)):
	data2=readPNSD_Mrg.readData(readThis=i,npsdStr=readPNSD_Mrg.npsdStr[i])
	data1[readPNSD_Mrg.npsdStr[i]]=data2[readPNSD_Mrg.npsdStr[i]].copy()

cdnc1=calc_average(cloud_times['Exp006'],\
	data1['MeteoCPC-Exp006'],data1['MergedOPC-Exp006'])
	
cdnc2=calc_average(cloud_times['Exp007'],\
	data1['MeteoCPC-Exp007'],data1['MergedOPC-Exp007'])
	
cdnc3=calc_average(cloud_times['Exp008'],\
	data1['MeteoCPC-Exp008'],data1['MergedOPC-Exp008'])

cdnc0=calc_average(cloud_times['Exp028'],\
	data1['MeteoCPC-Exp028'],data1['MergedOPC-Exp028'])
	
plt.ion()
plt.plot([0,600,2000],[cdnc0,cdnc1,cdnc2])
plt.ylim((0,800))