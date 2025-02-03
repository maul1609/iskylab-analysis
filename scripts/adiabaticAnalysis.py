import matplotlib.pyplot as plt
import numpy as np
import readMeteoCPC
import readOPC_Merged
import readModel
import svp
from scipy.interpolate import interp1d

timeCloud=[300, #2
	660,
	1080,
	380, #5
	540,
	690,
	580,
	510, #9
	500, #10
	660, #11
	530,  #12
	560,  #13
	480, 
	530,  #15
	490,  #16
	70, 
	360, #18
	500,  #19
	710,  #21
	]
Ra=287.
Rv=461.
eps1=Ra/Rv
def doPlot(readThisMet,readThisOPC,metStr,opcStr,title1,ax):
	dataMet=readMeteoCPC.readData(readThis=readThisMet,metStr=metStr)
	dataOPC=readOPC_Merged.readData(readThis=readThisOPC,opcStr=opcStr)
	
	
	# theta
	#theta=(dataMet[metStr]['Tgw mean']+273.15)*(1000.0/dataMet[metStr]['Pressure'])**(0.286)
	#plt.plot(dataMet[metStr]['Time'],theta)
	#plt.ylim((280,290))
	
	t0i=interp1d(dataMet[metStr]['Time'],\
	   dataMet[metStr]['Tgw mean'])
	t0=t0i(0)
	p0i=interp1d(dataMet[metStr]['Time'],\
	   dataMet[metStr]['Pressure'])
	p0=p0i(0)
	theta0=(t0+273.15)*(1000.0/p0)**0.286
	plot=4
	if plot==1:
		# plot actual temperature
		plt.plot(dataMet[metStr]['Time'],\
		   dataMet[metStr]['Tgw mean'])
		
		# plot calculated temperature - dry adiabatic
		plt.plot(dataMet[metStr]['Time'],\
		   theta0*(dataMet[metStr]['Pressure']/1000.0)**0.286-273.15)
	elif plot==2:
		#plt.plot(dataOPC[opcStr]['Time'],dataOPC[opcStr]['lwc'])
		# calculate the svmr at timeCloud
		tC=t0i(timeCloud[readThisMet])+273.15
		pC=np.array(p0i(timeCloud[readThisMet])*100.)
		eC=np.array(svp.svp([tC],'buck2','liq'))
		svmrC=eps1*eC/(pC-eC)
		svp1=np.array(svp.svp(dataMet[metStr]['Tgw mean']+273.15,'buck2','liq'))
		svmr=eps1*svp1/(dataMet[metStr]['Pressure']*100.-svp1)
		lwmr=np.maximum(svmrC-svmr,0.0)
		rhoa=dataMet[metStr]['Pressure']*100./(dataMet[metStr]['Tgw mean']+273.15)/Ra
		lwc=lwmr*rhoa
		dlwc=np.diff(lwc) # essentially dlwc / dt
		
		lwc_calc=np.zeros(len(lwc))
		for i in range(len(lwc)-1):
			lwc_calc[i+1]=dlwc[i]+lwc_calc[i]-1./50.*lwc_calc[i]
		plt.plot(dataOPC[opcStr]['Time'],lwc_calc*1000.)
		plt.plot(dataOPC[opcStr]['Time'],lwc*1000.)
		plt.plot(dataOPC[opcStr]['Time'],dataOPC[opcStr]['lwc'])
		plt.ylabel('LWC (g m$^{-3}$)')
		plt.ylim((0,2))
		#plt.plot(dataOPC[opcStr]['Time'],dataOPC[opcStr]['lwc']/(lwc*1000.))
		#plt.ylim((1e-3,1))
		#plt.ylabel('Adiabatic fraction')
		#plt.yscale('log')
		plt.grid()
	elif plot==3:
		svp1=np.array(svp.svp(dataMet[metStr]['Tgw mean']+273.15,'buck2','liq'))
		svmr=eps1*svp1/(dataMet[metStr]['Pressure']*100.-svp1)
		vp1=np.array(svp.svp(dataMet[metStr]['TDew']+273.15,'buck2','liq'))
		vmr=eps1*vp1/(dataMet[metStr]['Pressure']*100.-vp1)
		
		plt.plot(dataOPC[opcStr]['Time'],vmr)
		rhoa=dataMet[metStr]['Pressure']*100./(dataMet[metStr]['Tgw mean']+273.15)/Ra
		plt.plot(dataOPC[opcStr]['Time'],vmr+dataOPC[opcStr]['lwc']/rhoa/1000.)
		plt.grid()
	elif plot==4:
		#temp=theta0*(dataMet[metStr]['Pressure']/1000.)**0.286
		svp1=np.array(svp.svp(dataMet[metStr]['Tgw mean']+273.15,'buck2','liq'))
		svmr=eps1*svp1/(dataMet[metStr]['Pressure']*100.-svp1)
		vp1=np.array(svp.svp(dataMet[metStr]['TDew']+273.15,'buck2','liq'))
		vmr=eps1*vp1/(dataMet[metStr]['Pressure']*100.-vp1)
		vp1=np.array(svp.svp(dataMet[metStr]['Tgw mean']+273.15,'buck2','liq'))
		vmr2=eps1*vp1/(dataMet[metStr]['Pressure']*100.-vp1)
		vp1=np.array(svp.svp(dataMet[metStr]['Tww mean']+273.15,'buck2','liq'))
		vmr3=eps1*vp1/(dataMet[metStr]['Pressure']*100.-vp1)
		
		rhoa=dataMet[metStr]['Pressure']*100./(dataMet[metStr]['Tgw mean']+273.15)/Ra
		plt.plot(dataMet[metStr]['Time'],vmr)
		plt.plot(dataMet[metStr]['Time'],vmr[0]*np.ones(len(vmr)))
		vmr=vmr+dataOPC[opcStr]['lwc']/rhoa/1000.
		plt.plot(dataMet[metStr]['Time'],vmr)
		plt.plot(dataMet[metStr]['Time'],vmr2)
		plt.plot(dataMet[metStr]['Time'],vmr3)
		plt.ylim((0,0.0065))
		plt.plot([timeCloud[readThisMet],timeCloud[readThisMet]],plt.ylim())
		if(readThisMet==3):
			plt.legend(['actual vmr','initial vmr','actual vmr+lwmr','svmr at Tgas','svmr at Tw'])
		plt.grid()

	plt.text(0.1,0.8,title1,transform=ax.transAxes)
	#plt.xlim((-50,1500))

if __name__=="__main__":
	fig=plt.figure()
	plt.ion()
	"""
	ax=plt.subplot(121)
	# experiment 5, AS
	doPlot(26,26,'MeteoCPC-Exp028','MergedOPC-Exp028','Experiment 028: AS (0.01 wt%)',ax)
	ax=plt.subplot(122)
	# experiment 6, AS + NaCl
	doPlot(27,27,'MeteoCPC-Exp029','MergedOPC-Exp029','Experiment 029: AS (0.01 wt%)',ax)
	"""
	ax=plt.subplot(421)
	# experiment 5, AS
	doPlot(3,3,'MeteoCPC-Exp005','MergedOPC-Exp005','Experiment 005: AS (0.01 wt%)',ax)
	ax=plt.subplot(422)
	# experiment 6, AS + NaCl
	doPlot(4,4,'MeteoCPC-Exp006','MergedOPC-Exp006','Experiment 006: AS (0.01 wt%) + NaCl (0.1 wt%)',ax)
	ax=plt.subplot(423)
	# exoeriment 11, AS
	doPlot(9,9,'MeteoCPC-Exp011','MergedOPC-Exp011','Experiment 011: AS (0.01 wt%)',ax)
	ax=plt.subplot(424)
	# experiment 11, AS
	doPlot(9,9,'MeteoCPC-Exp011','MergedOPC-Exp011','Experiment 011: AS (10.0 wt%)',ax)
	ax=plt.subplot(425)
	# experiment 12, NaCl
	doPlot(10,10,'MeteoCPC-Exp012','MergedOPC-Exp012','Experiment 012: NaCl (0.1 wt%)',ax)
	ax=plt.subplot(426)
	# experiment 13, NaCl
	doPlot(11,11,'MeteoCPC-Exp013','MergedOPC-Exp013','Experiment 013: NaCl (0.1 wt%)',ax)
	ax=plt.subplot(427)
	# experiment 16, AS
	doPlot(14,14,'MeteoCPC-Exp016','MergedOPC-Exp016','Experiment 016: AS (1.0 wt%)',ax)
	ax=plt.subplot(428)
	# experiment 16, AS
	doPlot(15,15,'MeteoCPC-Exp017','MergedOPC-Exp017','Experiment 017: AS (1.0 wt%)',ax)
	
	fig.tight_layout()
	#plt.subplots_adjust(wspace=0, hspace=0)	
	plt.show()

