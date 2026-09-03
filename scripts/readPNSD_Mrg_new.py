"""Read and fit iSKYLAB initial particle-number size distributions.

Measured component PSDs are represented by up to three lognormal submodes.
Fitted diameters are kept in micrometres; conversion to SI metres is performed
only by the BMM namelist generator.
"""

import numpy as np
import matplotlib.pyplot as plt
import iskylab_config as cfg

# Fit each measured component PSD with up to three lognormal submodes.
# Diameters read/fitted here are in micrometres; the batch namelist generator
# performs the single conversion to metres when writing BMM aerosol input.
readThis = 7

fileNamePNSD_Mrg=[ \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp002-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp003-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp004-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp005-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp006-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp007-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp008-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp009-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp010-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp011-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp012-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp013-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp014-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp015-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp016-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp017-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp018-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp019-3-InitialPNSD-Mrg.csv'), \
str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp020-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp021-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp022-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp023-3-InitialPNSD-Mrg.csv'), \
# 	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp024-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp025-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp026-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp027-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp028-3-InitialPNSD-Mrg.csv'), \
	str(cfg.DATA_ROOT / 'Datasets-V2/Initial-PNSD/iSKYLAB01-Exp029-3-InitialPNSD-Mrg.csv')]

npsdStr=['InitialPNSD-Exp002','InitialPNSD-Exp003','InitialPNSD-Exp004','InitialPNSD-Exp005',\
	'InitialPNSD-Exp006','InitialPNSD-Exp007','InitialPNSD-Exp008','InitialPNSD-Exp009',\
	'InitialPNSD-Exp010','InitialPNSD-Exp011','InitialPNSD-Exp012','InitialPNSD-Exp013',\
	'InitialPNSD-Exp014','InitialPNSD-Exp015','InitialPNSD-Exp016','InitialPNSD-Exp017',\
	'InitialPNSD-Exp018','InitialPNSD-Exp019',\
	'InitialPNSD-Exp020','InitialPNSD-Exp021',\
	'InitialPNSD-Exp022','InitialPNSD-Exp023',\
	#'InitialPNSD-Exp024',
	'InitialPNSD-Exp025',\
	'InitialPNSD-Exp026','InitialPNSD-Exp027','InitialPNSD-Exp028','InitialPNSD-Exp029']

# note, got to include PSD parameters for SDSA01, SDTP02 and ATD03
whichPSD=[[1],[1],[5],[1],\
	[1,0],[1,0],[1,3],[1,3],[1,5],\
	[1],[0],[0],[0],[1],\
	[5],[5],[1,0],[1,0],\
	[-1],[-1],[-1,0],[-1,0],[-2],[-2,0],[-2,0],\
	[1],[1]]
targetConc=[[500],[5000],[5000],[3000],\
	[3000,600],[3000,2000],[3000,2000],[3000,4000],[3000,4000],\
	[5000],[2000],[600],[2000],[3000],[3000],[3000],[3000,600],[3000,2000],\
	[3000],[1000],[3000,2000],[3000,2000],[1000],[1000,4000],[1000,4000],\
	[3000],[3000]]
kappa=[[0.61],[0.61],[0.61],[0.61],\
	[0.61,1.28],[0.61,1.28],[0.61,1.28],[0.61,1.28],[0.61,1.28],\
	[0.61],[1.28],[1.28],[1.28],[0.067],[0.61],[0.61],[0.067,1.28],[0.067,1.28],\
	[0.0],[0.0],[0.0,1.28],[0.0,1.28],[0.0],[0.0,1.28],[0.0,1.28],\
	[0.61],[0.61]]
	
# kappa=[[0.53],[0.53],[0.53],[0.53],\
# 	[0.53,1.12],[0.53,1.12],[0.53,1.12],[0.53,1.12],[0.53,1.12],\
# 	[0.53],[1.12],[1.12],[1.12],[0.067],[0.53],[0.53],[0.067,1.12],[0.067,1.12],\
# 	[-1],[-1,1.12],[-1,1.12],[-2],[-2,1.12],[-2,1.12],\
# 	[0.53],[0.53]]

comp=[['AS'],['AS'],['AS'],['AS'],\
	['AS','NaCl'],['AS','NaCl'],['AS','NaCl'],['AS','NaCl'],['AS','NaCl'],\
	['AS'],['NaCl'],['NaCl'],['NaCl'],['FA'],['AS'],['AS'],['FA','NaCl'],['FA','NaCl'],\
	['SDSA01'],['SDSA01'],['SDSA01','NaCl'],['SDSA01','NaCl'],['ATD03'],['ATD03','NaCl'],['ATD03','NaCl'],\
	['AS'],['AS']]

density=[[1770],[1770],[1770],[1770],\
	[1770,2160],[1770,2160],[1770,2160],[1770,2160],[1770,2160],\
	[1770],[2160],[2160],[2160],[1500],[1770],[1770],[1500,2160],[1500,2160],\
	[4000],[4000],[4000,2160],[4000,2160],[4000],[4000,2160],[4000,2160],\
	[1770],[1770]]

N=[[0.49,0.38],[0.18,0.74],[0.16,0.91],[0.2,1.06],[0.6,1.37],[0.56,1.18]]
lnsig=[[0.25,0.84],[0.19,0.45],[0.19,0.43],[0.23,0.47],[0.49,0.76],[0.46,0.68]]
Dm=[[0.247,0.205],[0.122,0.14],[0.084,0.115],[0.061,0.102],[0.038,0.08],[0.029,0.053]]


# -----------------------------------------------------------------------------
# BMM-oriented PSD fitting
# -----------------------------------------------------------------------------
#
# BMM uses the fitted PSD for quantities other than total aerosol number:
#   * dry surface area (~D^2) controls IASD/INAS active-site abundance;
#   * dry volume (~D^3) controls aerosol mass;
#   * N(D >= 0.5 um) is important to DeMott-style diagnostics/prognostics.
#
# An unweighted least-squares fit to dN/dlnD is dominated by the number peak.
# The fitter below therefore combines a robust log-spectrum residual with
# explicit penalties on BMM-relevant moments.
PNSD_FIT_COARSE_THRESHOLD_UM = 0.5
PNSD_FIT_WEIGHTS = {
    "shape": 1.0,
    "number": 0.75,
    "surface_area": 3.0,
    "volume": 1.5,
    "coarse_number": 3.0,
}
PNSD_FIT_WARN_REL_ERROR = 0.10

# Broader than the historical bounds so coarse dust is not artificially forced
# to Dm <= 0.4 um or ln(sigma) <= 0.65.
PNSD_FIT_DMIN_UM = 0.005
PNSD_FIT_DMAX_UM = 2.0
PNSD_FIT_LNSIG_MIN = 0.08
PNSD_FIT_LNSIG_MAX = 1.20


def _single_lognormal_dndlnD(x, number, lnsig, dm):
    x = np.asarray(x, dtype=float)
    number = max(float(number), 0.0)
    lnsig = max(float(lnsig), np.finfo(float).tiny)
    dm = max(float(dm), np.finfo(float).tiny)
    return (
        number / (np.sqrt(2.0*np.pi)*lnsig)
        * np.exp(-(np.log(x/dm)**2.0)/(2.0*lnsig**2.0))
    )


def lognormal_func2(x, a,b,c,d,e,f,g,h,i):
    """Historical public 3-mode function retained for compatibility."""
    return (
        _single_lognormal_dndlnD(x,a,b,c)
        + _single_lognormal_dndlnD(x,d,e,f)
        + _single_lognormal_dndlnD(x,g,h,i)
    )


def _log_edges_from_centres(d):
    d=np.asarray(d,dtype=float)
    if d.ndim != 1 or d.size < 2:
        raise ValueError("PSD diameter must be a 1-D array with at least 2 bins")
    if np.any(~np.isfinite(d)) or np.any(d <= 0.0) or np.any(np.diff(d) <= 0.0):
        raise ValueError("PSD diameter bins must be finite, positive and increasing")
    ld=np.log(d)
    e=np.empty(d.size+1)
    e[1:-1]=0.5*(ld[:-1]+ld[1:])
    e[0]=ld[0]-0.5*(ld[1]-ld[0])
    e[-1]=ld[-1]+0.5*(ld[-1]-ld[-2])
    return e


def _measured_moments(d,y,threshold=PNSD_FIT_COARSE_THRESHOLD_UM):
    """Measured raw diameter moments using constant dN/dlnD in each log bin."""
    d=np.asarray(d,dtype=float)
    y=np.asarray(y,dtype=float)
    e=_log_edges_from_centres(d)
    dl=np.diff(e)

    out={}
    out["number"]=float(np.sum(y*dl))
    out["surface_area_proxy"]=float(np.sum(y*d**2*dl))
    out["volume_proxy"]=float(np.sum(y*d**3*dl))

    lthr=np.log(float(threshold))
    overlap=np.maximum(0.0,e[1:]-np.maximum(e[:-1],lthr))
    out["coarse_number"]=float(np.sum(y*overlap))
    return out


def _fitted_moments(number,lnsig,dm,threshold=PNSD_FIT_COARSE_THRESHOLD_UM):
    """Analytic number, D^2, D^3 and coarse-number moments of fitted modes."""
    from scipy.special import erfc

    number=np.asarray(number,dtype=float)
    lnsig=np.asarray(lnsig,dtype=float)
    dm=np.asarray(dm,dtype=float)

    out={}
    out["number"]=float(np.sum(number))
    out["surface_area_proxy"]=float(
        np.sum(number*dm**2*np.exp(2.0*lnsig**2))
    )
    out["volume_proxy"]=float(
        np.sum(number*dm**3*np.exp(4.5*lnsig**2))
    )
    z=np.log(float(threshold)/dm)/(np.sqrt(2.0)*lnsig)
    out["coarse_number"]=float(np.sum(0.5*number*erfc(z)))
    return out


def _fit_curve(d,number,lnsig,dm):
    y=np.zeros_like(np.asarray(d,dtype=float))
    for n,s,m in zip(number,lnsig,dm):
        y += _single_lognormal_dndlnD(d,n,s,m)
    return y


def _safe_log_ratio(a,b):
    tiny=1.0e-30
    return np.log(max(float(a),tiny)/max(float(b),tiny))


def _initial_guesses(d,y):
    """Independent deterministic multi-start guesses for one component."""
    e=_log_edges_from_centres(d)
    weights=np.maximum(y,0.0)*np.diff(e)
    total=float(np.sum(weights))
    if total <= 0.0:
        total=1.0

    cdf=np.cumsum(weights)
    if cdf[-1] > 0.0:
        cdf=cdf/cdf[-1]

    starts=[]
    for qs in ([0.15,0.50,0.85],[0.08,0.40,0.85],[0.25,0.65,0.93]):
        dm=[]
        for q in qs:
            k=int(np.searchsorted(cdf,q,side="left"))
            k=min(max(k,0),len(d)-1)
            dm.append(d[k])
        starts.append((
            np.full(3,total/3.0),
            np.full(3,0.35),
            np.asarray(dm,dtype=float),
        ))

    # Broad generic aerosol/dust starts.
    for dm,sig,frac in (
        ([0.03,0.10,0.35],[0.30,0.35,0.45],[0.45,0.40,0.15]),
        ([0.05,0.20,0.65],[0.35,0.45,0.55],[0.40,0.40,0.20]),
        ([0.08,0.30,0.90],[0.30,0.50,0.60],[0.35,0.45,0.20]),
    ):
        starts.append((
            total*np.asarray(frac),
            np.asarray(sig,dtype=float),
            np.asarray(dm,dtype=float),
        ))
    return starts


def _pack(number,lnsig,dm):
    return np.column_stack((number,lnsig,dm)).reshape(-1)


def _unpack(p):
    q=np.asarray(p,dtype=float).reshape(3,3)
    return q[:,0].copy(),q[:,1].copy(),q[:,2].copy()


def _fit_quality(d,y,number,lnsig,dm):
    fitted=_fit_curve(d,number,lnsig,dm)
    positive=y[y > 0.0]
    floor=max(np.percentile(positive,10.0)*0.1,1.0e-6)
    log_rmse=float(np.sqrt(np.mean(
        (np.log10(fitted+floor)-np.log10(y+floor))**2
    )))

    measured=_measured_moments(d,y)
    model=_fitted_moments(number,lnsig,dm)
    rel={}
    for name,val in measured.items():
        rel[name]=(model[name]-val)/val if val > 0.0 else np.nan

    return {
        "log10_rmse":log_rmse,
        "measured":measured,
        "fitted":model,
        "relative_error":rel,
    }


def _fit_component_bmm(diameter_um,observed):
    """Three-lognormal fit optimised for the quantities BMM actually uses."""
    from scipy.optimize import least_squares

    d=np.asarray(diameter_um,dtype=float)
    y=np.asarray(observed,dtype=float)

    valid=np.isfinite(d) & np.isfinite(y) & (d > 0.0) & (y >= 0.0)
    d=d[valid]
    y=y[valid]
    if d.size < 6:
        raise ValueError("Too few valid PSD bins for a 3-mode fit")
    if not np.any(y > 0.0):
        raise ValueError("PSD contains no positive concentrations")

    measured=_measured_moments(d,y)

    positive=y[y > 0.0]
    floor=max(np.percentile(positive,10.0)*0.1,1.0e-6)

    nscale=max(measured["number"],1.0)
    nmax=max(20.0*nscale,1.0e5)

    lower=np.tile(
        [0.0,PNSD_FIT_LNSIG_MIN,PNSD_FIT_DMIN_UM],3
    )
    upper=np.tile(
        [nmax,PNSD_FIT_LNSIG_MAX,PNSD_FIT_DMAX_UM],3
    )

    shape_scale=PNSD_FIT_WEIGHTS["shape"]/np.sqrt(d.size)

    def residual(p):
        number,lnsig,dm=_unpack(p)
        model_y=_fit_curve(d,number,lnsig,dm)
        r=list(
            shape_scale*(np.log(model_y+floor)-np.log(y+floor))
        )

        fm=_fitted_moments(number,lnsig,dm)
        r.append(
            PNSD_FIT_WEIGHTS["number"]*
            _safe_log_ratio(fm["number"],measured["number"])
        )
        r.append(
            PNSD_FIT_WEIGHTS["surface_area"]*
            _safe_log_ratio(
                fm["surface_area_proxy"],measured["surface_area_proxy"]
            )
        )
        r.append(
            PNSD_FIT_WEIGHTS["volume"]*
            _safe_log_ratio(fm["volume_proxy"],measured["volume_proxy"])
        )
        if measured["coarse_number"] > 0.0:
            r.append(
                PNSD_FIT_WEIGHTS["coarse_number"]*
                _safe_log_ratio(
                    fm["coarse_number"],measured["coarse_number"]
                )
            )
        return np.asarray(r)

    best=None
    for nums,sigs,dms in _initial_guesses(d,y):
        x0=_pack(
            np.clip(nums,1.0e-10,0.9*nmax),
            np.clip(sigs,PNSD_FIT_LNSIG_MIN*1.01,PNSD_FIT_LNSIG_MAX*0.99),
            np.clip(dms,PNSD_FIT_DMIN_UM*1.01,PNSD_FIT_DMAX_UM*0.99),
        )
        result=least_squares(
            residual,x0,bounds=(lower,upper),
            method="trf",loss="soft_l1",f_scale=0.5,
            x_scale="jac",max_nfev=5000,
        )
        if best is None or result.cost < best.cost:
            best=result

    if best is None or not best.success:
        raise RuntimeError("BMM-oriented PSD fit failed")

    number,lnsig,dm=_unpack(best.x)
    order=np.argsort(dm)
    number=number[order]
    lnsig=lnsig[order]
    dm=dm[order]

    metrics=_fit_quality(d,y,number,lnsig,dm)
    metrics["optimiser_cost"]=float(best.cost)
    metrics["optimiser_nfev"]=int(best.nfev)

    # Flag active modes very close to a physical bound.
    bound_hits=[]
    tol=2.0e-3
    for k,(n,s,m) in enumerate(zip(number,lnsig,dm),start=1):
        if n <= 1.0e-6*nscale:
            continue
        hits=[]
        if abs(s-PNSD_FIT_LNSIG_MIN) < tol: hits.append("lnsig_min")
        if abs(s-PNSD_FIT_LNSIG_MAX) < tol: hits.append("lnsig_max")
        if abs(m-PNSD_FIT_DMIN_UM) < tol: hits.append("Dm_min")
        if abs(m-PNSD_FIT_DMAX_UM) < tol: hits.append("Dm_max")
        if hits:
            bound_hits.append((k,hits))
    metrics["bound_hits"]=bound_hits

    return number,lnsig,dm,metrics


def _print_fit_quality(exp_name,component_key,metrics):
    rel=metrics["relative_error"]

    def pct(name):
        x=rel.get(name,np.nan)
        return "n/a" if not np.isfinite(x) else f"{100.0*x:+.2f}%"

    print(
        f"{exp_name} {component_key}: "
        f"log10-RMSE={metrics['log10_rmse']:.3f}; "
        f"N={pct('number')}; "
        f"area={pct('surface_area_proxy')}; "
        f"volume={pct('volume_proxy')}; "
        f"N(D>={PNSD_FIT_COARSE_THRESHOLD_UM:g}um)="
        f"{pct('coarse_number')}"
    )

    important=("surface_area_proxy","volume_proxy","coarse_number")
    bad=[
        name for name in important
        if np.isfinite(rel.get(name,np.nan))
        and abs(rel[name]) > PNSD_FIT_WARN_REL_ERROR
    ]
    if bad:
        print(
            "  WARNING: >"
            f"{100.0*PNSD_FIT_WARN_REL_ERROR:.0f}% BMM-moment error in "
            + ", ".join(bad)
        )
    if metrics["bound_hits"]:
        print("  WARNING: fitted mode near bound:",metrics["bound_hits"])


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

	keys1=data1[npsdStr].keys()
	off1=2
	num1=int(float(len(list(keys1)[2:]))/2)
	keyList=list(keys1)

	# change the key name if only one
	if(num1==1):
		temp=data1[npsdStr][keyList[off1+num1]]
		del data1[npsdStr][keyList[off1+num1]]
		keyList[off1+num1]=keyList[off1+num1].replace('Total',comp[readThis][0])
		data1[npsdStr][keyList[off1+num1]]=temp.copy()


	data1[npsdStr]['num1']=num1
	data1[npsdStr]['keyList']=keyList.copy()

	for j in range(num1):
		component_key=keyList[off1+num1+j]
		observed=data1[npsdStr][component_key]

		N2,lnsig2,dm2,metrics=_fit_component_bmm(
			data1[npsdStr]['Dve'],observed
		)

		data1[npsdStr]['Nfit_' + component_key]=N2.tolist()
		data1[npsdStr]['lnsigfit_' + component_key]=lnsig2.tolist()
		data1[npsdStr]['dfit_' + component_key]=dm2.tolist()

		# Store diagnostics next to the fitted parameters so downstream scripts
		# can inspect fit quality without refitting the PSD.
		data1[npsdStr]['fit_metrics_' + component_key]=metrics
		_print_fit_quality(npsdStr,component_key,metrics)

	
	return data1

if __name__ == "__main__":
	doAnalysis = True
	off1=2
	
	d=np.logspace(-2,np.log10(2),100)
	
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
		keyList=data1[npsdStr[i]]['keyList']
		num1=data1[npsdStr[i]]['num1']
		for j in range(num1):
			plt.plot(data1[npsdStr[i]]['Dve'], \
				data1[npsdStr[i]][keyList[off1+num1+j]])
		
		for j in range(num1):
			""" 
				do the fit
			"""			
			N2=data1[npsdStr[i]]['Nfit_' + keyList[off1+num1+j]]
			lnsig2=data1[npsdStr[i]]['lnsigfit_' + keyList[off1+num1+j]]
			dm2=data1[npsdStr[i]]['dfit_' + keyList[off1+num1+j]]
			dNdlogD=np.zeros(len(d))
			
			for k in range(len(dm2)):
				dNdlogD=dNdlogD+ \
					N2[k]/(np.sqrt(2.0*np.pi)*lnsig2[k])* \
					np.exp(-(np.log(d/dm2[k])**2.0)/(2*lnsig2[k]**2))
			plt.plot(d,dNdlogD,lw=0.5,color='k')
			
		plt.legend(keyList[off1+num1:])
		plt.title(npsdStr[i])

		plt.xscale('log')
		plt.xlim((0.01,2))
