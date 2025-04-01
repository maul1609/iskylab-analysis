import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# derive fits for each individual aeroosl PSD
# calculate ratios for each case
# for each individual fit recalculate based on slip correction
# now you have initial conditions
readThis = 7

fileNamePNSD_Mrg=[ \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp002-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp003-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp004-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp005-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp006-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp007-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp008-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp009-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp010-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp011-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp012-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp013-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp014-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp015-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp016-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp017-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp018-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp019-3-InitialPNSD-Mrg.csv', \
# 	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp020-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp021-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp022-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp023-3-InitialPNSD-Mrg.csv', \
# 	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp024-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp025-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp026-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp027-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp028-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/Initial-PNSD/iSKYLAB01-Exp029-3-InitialPNSD-Mrg.csv']

npsdStr=['InitialPNSD-Exp002','InitialPNSD-Exp003','InitialPNSD-Exp004','InitialPNSD-Exp005',\
	'InitialPNSD-Exp006','InitialPNSD-Exp007','InitialPNSD-Exp008','InitialPNSD-Exp009',\
	'InitialPNSD-Exp010','InitialPNSD-Exp011','InitialPNSD-Exp012','InitialPNSD-Exp013',\
	'InitialPNSD-Exp014','InitialPNSD-Exp015','InitialPNSD-Exp016','InitialPNSD-Exp017',\
	'InitialPNSD-Exp018','InitialPNSD-Exp019',\
	#'InitialPNSD-Exp020', 
	'InitialPNSD-Exp021',\
	'InitialPNSD-Exp022','InitialPNSD-Exp023',\
	#'InitialPNSD-Exp024',
	'InitialPNSD-Exp025',\
	'InitialPNSD-Exp026','InitialPNSD-Exp027','InitialPNSD-Exp028','InitialPNSD-Exp029']

# note, got to include PSD parameters for SDSA01, SDTP02 and ATD03
whichPSD=[[1],[1],[5],[1],\
	[1,0],[1,0],[1,3],[1,3],[1,5],\
	[1],[0],[0],[0],[1],\
	[5],[5],[1,0],[1,0],\
	[-1],[-1,0],[-1,0],[-2],[-2,0],[-2,0],\
	[1],[1]]
targetConc=[[500],[5000],[5000],[3000],\
	[3000,600],[3000,2000],[3000,2000],[3000,4000],[3000,4000],\
	[5000],[2000],[600],[2000],[3000],[3000],[3000],[3000,600],[3000,2000],\
	[1000],[3000,2000],[3000,2000],1000,[1000,4000],[1000,10000],\
	[3000],[3000]]
kappa=[[0.61],[0.61],[0.61],[0.61],\
	[0.61,1.28],[0.61,1.28],[0.61,1.28],[0.61,1.28],[0.61,1.28],\
	[0.61],[1.28],[1.28],[1.28],[0.067],[0.61],[0.61],[0.067,1.28],[0.067,1.28],\
	[-1],[-1,1.28],[-1,1.28],[-2],[-2,1.28],[-2,1.28],\
	[0.61],[0.61]]
density=[[1770],[1770],[1770],[1770],\
	[1770,2160],[1770,2160],[1770,2160],[1770,2160],[1770,2160],\
	[1770],[2160],[2160],[2160],[1500],[1770],[1770],[1500,2160],[1500,2160],\
	[4000],[4000,2160],[4000,2160],[4000],[4000,2160],[4000,2160],\
	[1770],[1770]]

N=[[0.49,0.38],[0.18,0.74],[0.16,0.91],[0.2,1.06],[0.6,1.37],[0.56,1.18]]
lnsig=[[0.25,0.84],[0.19,0.45],[0.19,0.43],[0.23,0.47],[0.49,0.76],[0.46,0.68]]
Dm=[[0.247,0.205],[0.122,0.14],[0.084,0.115],[0.061,0.102],[0.038,0.08],[0.029,0.053]]


# Define the nonlinear function
def lognormal_func2(x, a,b,c,d,e,f,g,h,i):
	
	dNdlogD= \
		a/(np.sqrt(2.0*np.pi)*b)* \
		np.exp(-(np.log(x/c)**2.0)/(2*b**2))
	dNdlogD=dNdlogD+ \
		d/(np.sqrt(2.0*np.pi)*e)* \
		np.exp(-(np.log(x/f)**2.0)/(2*e**2))
	dNdlogD=dNdlogD+ \
		g/(np.sqrt(2.0*np.pi)*h)* \
		np.exp(-(np.log(x/i)**2.0)/(2*h**2))
# 	print(a,b,c,d,e,f)
	return dNdlogD



def readData(readThis = 3,npsdStr="InitialPNSD-Exp005"):
	fp = open(fileNamePNSD_Mrg[readThis],'r')
	str1=fp.readlines()
	fp.close()
	
	"""
		diameters
	"""
	dtemp=str1[0].split(',')
	Dve=[float(val) for val in dtemp[1:-1]]
	dtemp=str1[1].split(',')
	dlogDve=[float(val) for val in dtemp[1:]]
	
	
	data1=dict()
	data1 = {npsdStr : \
		{"Dve" : np.array(Dve), "dlogDve" : np.array(dlogDve)}}
	for i in range(len(str1[2:])):
		dtemp=str1[2+i].split(',')
		name1=dtemp[0].replace('(','') \
			.replace(')','').replace('#','') \
			.replace('Dve','Dve_').replace('dN/dlog','_dN/dlog') \
			.replace(' ','').replace('/','')
		var1=[float(val)/np.log(10.0) for val in dtemp[1:-1]]
		data1[npsdStr][name1]=np.array(var1)

	return data1

if __name__ == "__main__":
	doAnalysis = True
	
	
	if 'data1' in locals():
		pass
	else:
		data1=dict()

	for i in range(len(npsdStr)):
		data2=readData(readThis=i,npsdStr=npsdStr[i])
		data1[npsdStr[i]]=data2[npsdStr[i]].copy()
		"""	
	if doAnalysis:
		"""
		plt.ion()
		plt.figure()
		keys1=data1[npsdStr[i]].keys()
		off1=2
		num1=int(float(len(list(keys1)[2:]))/2)
		keyList=list(keys1)
		d=np.logspace(-2,np.log10(2),100)
		dm2=[0.26,0.05,0.2]
		lnsig2=[0.2,0.2,0.3]
		N2=[3000*0.6, \
			3000*0.4,1000.]	
		
		
		for j in range(num1):
			plt.plot(data1[npsdStr[i]]['Dve'], \
				data1[npsdStr[i]][keyList[off1+num1+j]])
		
		
		for j in range(num1):
			""" 
				do the fit
			"""
# 			ind,=np.where(data1[npsdStr[readThis]][keyList[off1+num1+j]]>0.0)
			ind=np.mgrid[0:len(data1[npsdStr[i]][keyList[off1+num1+j]])]
			popt, pcov = curve_fit(lognormal_func2, data1[npsdStr[i]]['Dve'][ind],\
				data1[npsdStr[i]][keyList[off1+num1+j]][ind], \
				p0=[N2[0], lnsig2[0], dm2[0], N2[1], lnsig2[1],dm2[1], N2[2], lnsig2[2],dm2[2]],\
				bounds=([0.1,0.1,0.03,0.1,0.1,0.03,0.1,0.1,0.03],\
					[5000,0.65,0.4,5000,0.65,0.4,5000,0.65,0.4]), \
				method='trf') 
			dNdlogD=np.zeros(len(d))
			N2=[popt[0],popt[3],popt[6]]
			lnsig2=[popt[1],popt[4],popt[7]]
			dm2=[popt[2],popt[5],popt[8]]
			
			data1[npsdStr[i]]['Nfit_' + keyList[off1+num1+j]]=N2.copy()
			data1[npsdStr[i]]['lnsigfit_' + keyList[off1+num1+j]]=lnsig2.copy()
			data1[npsdStr[i]]['dfit_' + keyList[off1+num1+j]]=dm2.copy()
			
			for k in range(len(dm2)):
				dNdlogD=dNdlogD+ \
					N2[k]/(np.sqrt(2.0*np.pi)*lnsig2[k])* \
					np.exp(-(np.log(d/dm2[k])**2.0)/(2*lnsig2[k]**2))
			plt.plot(d,dNdlogD,lw=0.5,color='k')
			
		plt.legend(keyList[off1+num1:])
		plt.title(npsdStr[i])

		plt.xscale('log')
		plt.xlim((0.01,2))
