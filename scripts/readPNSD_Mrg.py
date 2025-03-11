import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# derive fits for each individual aeroosl PSD
# calculate ratios for each case
# for each individual fit recalculate based on slip correction
# now you have initial conditions
readThis = 9

fileNamePNSD_Mrg=[ \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp002-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp003-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp004-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp005-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp006-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp007-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp008-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp009-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp010-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp011-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp012-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp013-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp014-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp015-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp016-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp017-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp018-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp019-3-InitialPNSD-Mrg.csv', \
# 	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp020-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp021-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp022-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp023-3-InitialPNSD-Mrg.csv', \
# 	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp024-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp025-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp026-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp027-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp028-3-InitialPNSD-Mrg.csv', \
	'../iSKYLAB-data/Datasets/iSKYLAB01-Exp029-3-InitialPNSD-Mrg.csv']

# note, got to include PSD parameters for SDSA01, SDTP02 and ATD03
whichPSD=[[1],[1],[5],[1],\
	[1,0],[1,0],[1,3],[1,3],[1,5],\
	[1],[0],[0],[0],[1],\
	[1],[1],[1,0],[1,0]]
N=[[0.49,0.38],[0.18,0.74],[0.16,0.91],[0.2,1.06],[0.6,1.37],[0.56,1.18]]
lnsig=[[0.25,0.84],[0.19,0.45],[0.19,0.43],[0.23,0.47],[0.49,0.76],[0.46,0.68]]
Dm=[[0.247,0.205],[0.122,0.14],[0.084,0.115],[0.061,0.102],[0.038,0.08],[0.029,0.053]]

Dm[0]=[0.22,0.17]
lnsig[1]=[0.19,0.40]
Dm[1]=[0.122,0.13]

# Define the nonlinear function
def lognormal_func2(x, a,b,c,d,e,f):
	
	dNdlogD= \
		a/(np.sqrt(2.0*np.pi)*b)* \
		np.exp(-(np.log(x/c)**2.0)/(2*b**2))
	dNdlogD=dNdlogD+ \
		d/(np.sqrt(2.0*np.pi)*e)* \
		np.exp(-(np.log(x/f)**2.0)/(2*e**2))
	return dNdlogD

# Define the nonlinear function
def lognormal_func(x, *args):
	print(args)
	args1=[i for i in args] 
	print(args1)
	if(len(np.shape(args1))>1):
		args1=args1[0]
	dNdlogD=np.zeros_like(x)
	for j in range(len(whichPSD[readThis])):
		NThis=N[whichPSD[readThis][j]]
		lnsigThis=lnsig[whichPSD[readThis][j]]
		DmThis=Dm[whichPSD[readThis][j]]
		for i in range(len(NThis)):
			dNdlogD=dNdlogD + \
				args1[j] * \
				NThis[i] / \
				(np.sqrt(2.0*np.pi)*lnsigThis[i])* \
				np.exp(-(np.log(x/DmThis[i])**2.0) / \
				(2*lnsigThis[i]**2))
		
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
	dtemp=str1[2].split(',')
	dN_dlogDve_scc=[float(val) for val in dtemp[1:-1]]
	dtemp=str1[3].split(',')
	dN_dlogDve_cc=[float(val) for val in dtemp[1:-1]]
	dN_dlogDve_cc=dN_dlogDve_cc/np.log(10.0)

	data1=dict()
	data1 = {npsdStr : \
		{"Dve" : np.array(Dve), "dlogDve" : np.array(dlogDve), \
		"dN_dlogDve_scc" : np.array(dN_dlogDve_scc), \
		"dN_dlogDve_cc": np.array(dN_dlogDve_cc)}}
	return data1

if __name__ == "__main__":
	doAnalysis = True
	
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
	
	if 'data1' in locals():
		pass
	else:
		data1=dict()

	for i in range(len(npsdStr)):
		data2=readData(readThis=i,npsdStr=npsdStr[i])
		data1[npsdStr[i]]=data2[npsdStr[i]].copy()
	
	if doAnalysis:
		plt.ion()
		plt.figure()
		plt.plot(data1[npsdStr[readThis]]['Dve'],data1[npsdStr[readThis]]['dN_dlogDve_cc'])
		plt.xscale('log')
		plt.xlim((0.01,2))
		type1=2
		
		if type1 == 2:
			d=np.logspace(-2,0,100)
			dm2=[0.26,0.22]
			lnsig2=[0.2,0.6]
			N2=[3000*0.6*np.sqrt(2.0*np.pi)*lnsig2[0], \
				3000*0.4*np.sqrt(2.0*np.pi)*lnsig2[1]]	
			popt, pcov = curve_fit(lognormal_func2, data1[npsdStr[readThis]]['Dve'],\
				data1[npsdStr[readThis]]['dN_dlogDve_cc'], \
				p0=[N2[0], lnsig2[0], dm2[0], N2[1], lnsig2[1],dm2[1]],\
				method='trf') 
			dNdlogD=np.zeros(len(d))
			N2=[popt[0],popt[3]]
			lnsig2=[popt[1],popt[4]]
			dm2=[popt[2],popt[5]]
			for i in range(len(dm2)):
				dNdlogD=dNdlogD+ \
					N2[i]/(np.sqrt(2.0*np.pi)*lnsig2[i])* \
					np.exp(-(np.log(d/dm2[i])**2.0)/(2*lnsig2[i]**2))
			plt.plot(d,dNdlogD)
			
		else:
			popt, pcov = curve_fit(lognormal_func,data1[npsdStr[readThis]]['Dve'], \
				data1[npsdStr[readThis]]['dN_dlogDve_cc'], \
				p0=[3000.0,600], method='lm') 
			plt.plot(data1[npsdStr[readThis]]['Dve'], \
				lognormal_func(np.array(data1[npsdStr[readThis]]['Dve']),popt))
		
