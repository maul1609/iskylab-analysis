#!/bin/bash

###### WARM EXPERIMENTS #########
python3 dataAnalysis_new.py --experiment Exp005 --cloud-formation-time-min 8 \
	--bl-tau-s 80 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.35   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp006 --cloud-formation-time-min 12 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.09   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp007 --cloud-formation-time-min 12 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.1   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp008 --cloud-formation-time-min 12 \
	--bl-tau-s 10 --bl-evap-size-exp 1.2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.095   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp009 --cloud-formation-time-min 12 \
	--bl-tau-s 10 --bl-evap-size-exp 1.2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.095   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp010 --cloud-formation-time-min 11 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.095   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp011 --cloud-formation-time-min 13 \
	--bl-tau-s 10 --bl-evap-size-exp 1.5 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.1   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

# 12 seems like issues

python3 dataAnalysis_new.py --experiment Exp013 --cloud-formation-time-min 11 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.085   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp014 --cloud-formation-time-min 10.5 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.106   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

# FA is a bit strage - multiple activation
python3 dataAnalysis_new.py --experiment Exp015 --cloud-formation-time-min 10.75 \
	--bl-tau-s 10 --bl-evap-size-exp 1.2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.075   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion
	
python3 dataAnalysis_new.py --experiment Exp016 --cloud-formation-time-min 10.6 \
	--bl-tau-s 10 --bl-evap-size-exp 1.75 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.08   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

# exp 17 evidence that larger wall temperature difference makes a difference
python3 dataAnalysis_new.py --experiment Exp017 --cloud-formation-time-min 1.25 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.55   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp018 --cloud-formation-time-min 7 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.09   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

python3 dataAnalysis_new.py --experiment Exp019 --cloud-formation-time-min 11 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k -0.1   --wall-vapour-transfer-velocity-ms .5e-3 --plt-ion

# MIXED-PHASE
#python3 dataAnalysis_new.py --experiment Exp021 --cloud-formation-time-min 12. \
#	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
#	--bl-temp-offset-k 0.4    --plt-ion

# note, had to multiply Niemand by around 5.
python3 dataAnalysis_new.py --experiment Exp021  --initial-rh 0.70 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k 0.4   --iasd-xthresh 0.1   --plt-ion


# for the rest of the mixed-phase experiments I think I need to have the option to specify initial
# RH instead of cloud onset time. The reason is because I think a lot of them do not 100% saturate
# also an option that INAS particles do not need to be "activated" in order to freeze. They
# probably just need some liquid on them , but how much?

# MIXED-PHASE SENSITIVITY EXAMPLES
#
# If the measured trajectory never reaches water saturation, bypass the
# cloud-formation-time inference and prescribe RH at t=0 directly:
#
# python3 dataAnalysis_new.py --experiment Exp022 --initial-rh 0.95 \
#     --bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0
#
# To let IASD/INAS aerosol freeze before conventional activation, prescribe
# Xthresh.  This writes inas_wetting_mode=1 and inas_xthresh=Xthresh:
#
# python3 dataAnalysis_new.py --experiment Exp022 --initial-rh 0.95 \
#     --iasd-xthresh 0.04 \
#     --bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0

# note, had to multiply Niemand by around 5.
python3 dataAnalysis_new.py --experiment Exp022  --initial-rh 0.70 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k 0.4   --iasd-xthresh 0.1   --plt-ion

# not great had to multiply by 15
python3 dataAnalysis_new.py --experiment Exp023  --initial-rh 0.65 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k 0.4   --iasd-xthresh 0.1   --plt-ion
	
# ATD
python3 dataAnalysis_new.py --experiment Exp025  --initial-rh 0.70 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k 0.4   --iasd-xthresh 0.1   --plt-ion
	
python3 dataAnalysis_new.py --experiment Exp026  --initial-rh 0.70 \
	--bl-tau-s 10 --bl-evap-size-exp 2 --bl-wall-water-mode 0  \
	--bl-temp-offset-k 0.4   --iasd-xthresh 0.1   --plt-ion


