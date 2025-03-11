import matplotlib.pyplot as plt
import numpy as np
import readMeteoCPC
import readOPC_Merged
import readModel
import svp

def doPlot(readThisMet,readThisOPC,metStr,opcStr,modelFile,modelStr,title1):
	dataMet=readMeteoCPC.readData(readThis=readThisMet,metStr=metStr)
	dataOPC=readOPC_Merged.readData(readThis=readThisOPC,opcStr=opcStr)
	dataModel=readModel.readData(modelFile,modelStr)
	
	plt.ion()
	fig = plt.figure(figsize=(8,10))
	
	# pressure plot
	ax=plt.subplot(611)
	ax.plot(dataMet[metStr]['Time'],dataMet[metStr]['Pressure'])
	ax.plot(dataModel[modelStr]['time'],dataModel[modelStr]['p']/100.,'--')
	ax.set_yticks(np.mgrid[680:1100:40])
	ax.set_ylim((690,1000))
	ax.grid()
	ax.set_ylabel('Pressure (hPa)')
	ax.text(0.4,0.1,'(a) Pressure',transform=ax.transAxes)
	ax.legend(['AIDA','BMM'])
	

	# temperature plot
	ax=plt.subplot(612)
	ax.plot(dataMet[metStr]['Time'],dataMet[metStr]['Tgw mean'])
	ax.plot(dataModel[modelStr]['time'],dataModel[modelStr]['t']-273.15,'--')
	ax.set_yticks(np.mgrid[-10:10:2])
	ax.set_ylim((-12,6))
	ax.grid()
	ax.set_ylabel('Temperature ($^\circ$C)')
	ax.text(0.4,0.1,'(b) Temperature',transform=ax.transAxes)
	ax.legend(['AIDA','BMM'])
	xl=ax.get_xlim()

	# number concentration plot
	ax=plt.subplot(613)
	ax.plot(dataMet[metStr]['Time'],dataMet[metStr]['CPC_TotBot'])
	ax.plot(dataOPC[opcStr]['Time'],dataOPC[opcStr]['ndrop'])
	ax.plot(dataModel[modelStr]['time'],dataModel[modelStr]['ndrop']/1e6* \
		dataModel[modelStr]['rhoa'],'--')
	ax.set_yticks(np.mgrid[0:10000:500])
	ax.set_ylim((-500,4000))
	ax.set_xlim(xl)
	ax.grid()
	ax.set_ylabel('Numb. Conc. (cm$^{-3}$)')
	ax.text(0.4,0.9,'(c) Number concs.',transform=ax.transAxes)
	ax.legend(['AIDA CPC','AIDA Drops','BMM Drops'])


	# size plot
	ax=plt.subplot(614)
	ax.plot(dataOPC[opcStr]['Time'],dataOPC[opcStr]['Deff'])
	ax.plot(dataModel[modelStr]['time'],dataModel[modelStr]['deff']*1e6,'--')
	ax.set_yticks(np.mgrid[0:20:2])
	ax.set_ylim((-2,16))
	ax.set_xlim(xl)
	ax.grid()
	ax.set_ylabel('$D_{eff}$ ($\mu$m)')
	ax.text(0.4,0.9,'(d) Effective diam.',transform=ax.transAxes)
	ax.legend(['AIDA','BMM'])


	# wc plot
	ax=plt.subplot(615)
	ax.plot(dataOPC[opcStr]['Time'],dataOPC[opcStr]['lwc'])
	ax.plot(dataModel[modelStr]['time'],dataModel[modelStr]['ql']*1e3* \
		dataModel[modelStr]['rhoa'],'--')
	
	# add in MBW
	e=svp.svp(dataMet[metStr]['TDew']+273.15,'buck2','liq')
	esat=np.array(svp.svp(dataMet[metStr]['Tgw mean']+273.15,'buck2','liq'))*0.95
	p=dataMet[metStr]['Pressure']*100.0
	ax.plot(dataMet[metStr]['Time'], 0.622*(e-esat)/p*1000.0 )	
	
	ax.set_yticks(np.mgrid[0:2:0.2])
	ax.set_ylim((-0.2,2))
	ax.set_xlim(xl)
	ax.grid()
	ax.set_ylabel('LWC (g m$^{-3}$)')
	ax.text(0.4,0.9,'(e) Liquid water content',transform=ax.transAxes)
	ax.legend(['AIDA','BMM','MBW'])


	# rh plot
	ax=plt.subplot(616)
	e=svp.svp(dataMet[metStr]['TDew']+273.15,'buck2','liq')
	esat=svp.svp(dataMet[metStr]['Tgw mean']+273.15,'buck2','liq')
	p=dataMet[metStr]['Pressure']*100.0
	ax.plot(dataMet[metStr]['Time'], \
		e/(p-e) / (esat/(p-esat)) )
	ax.plot(dataModel[modelStr]['time'],dataModel[modelStr]['rh'],'--')
# 	plt.yticks(np.mgrid[0:20:2])
# 	plt.ylim((0,16))
	ax.set_xlim(xl)
	ax.grid()
	ax.set_ylabel('$RH_{liq}$')
	ax.set_xlabel('Time (s)')
	ax.text(0.4,0.1,'(f) Relative Humidity',transform=ax.transAxes)
	ax.legend(['AIDA','BMM'])



	plt.suptitle(title1)
	fig.tight_layout()
	plt.subplots_adjust(wspace=0, hspace=0)
	
	plt.show()

if __name__=="__main__":
	# experiment 5, AS
	doPlot(3,3,'MeteoCPC-Exp005','MergedOPC-Exp005', \
		'../iSKYLAB-data/modelOutput/outputExp005.nc','model-Exp005',\
		'Experiment 005: AS (0.01 wt%)')
	plt.savefig('/tmp/exp5.png')
	# experiment 6, AS + NaCl
	doPlot(4,4,'MeteoCPC-Exp006','MergedOPC-Exp006', \
		'../iSKYLAB-data/modelOutput/outputExp006a.nc','model-Exp006',\
		'Experiment 006: AS (0.01 wt%) + NaCl (0.1 wt%)')
	plt.savefig('/tmp/exp6.png')
	# exoeriment 11, AS
	doPlot(9,9,'MeteoCPC-Exp011','MergedOPC-Exp011', \
		'../iSKYLAB-data/modelOutput/outputExp011.nc','model-Exp011',\
		'Experiment 011: AS (0.01 wt%)')
	# experiment 11, AS
	doPlot(9,9,'MeteoCPC-Exp011','MergedOPC-Exp011', \
		'../iSKYLAB-data/modelOutput/outputExp011-2.nc','model-Exp011',\
		'Experiment 011: AS (10.0 wt%)')
	# experiment 11, AS with vapour sink
# 	doPlot(9,9,'MeteoCPC-Exp011','MergedOPC-Exp011', \
# 		'../modelOutput/outputExp011-v1e-3.nc','model-Exp011',\
# 		'Experiment 011: AS (0.01 wt% with vapour sink)')
	# experiment 12, NaCl
	doPlot(10,10,'MeteoCPC-Exp012','MergedOPC-Exp012', \
		'../iSKYLAB-data/modelOutput/outputExp012.nc','model-Exp012',\
		'Experiment 012: NaCl (0.1 wt%)')
	# experiment 13, NaCl
	doPlot(11,11,'MeteoCPC-Exp013','MergedOPC-Exp013', \
		'../iSKYLAB-data/modelOutput/outputExp013.nc','model-Exp013',\
		'Experiment 013: NaCl (0.1 wt%)')
	# experiment 16, AS
	doPlot(14,14,'MeteoCPC-Exp016','MergedOPC-Exp016', \
		'../iSKYLAB-data/modelOutput/outputExp016.nc','model-Exp016',\
		'Experiment 016: AS (1.0 wt%)')
	# experiment 16, AS
	doPlot(15,15,'MeteoCPC-Exp017','MergedOPC-Exp017', \
		'../iSKYLAB-data/modelOutput/outputExp017.nc','model-Exp017',\
		'Experiment 017: AS (1.0 wt%)')
