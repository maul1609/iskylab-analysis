Four categories of files are provided for each experiment.

    CampaignName-ExperimentNumber-0-Info.txt
        Basic information for experiment

    CampaignName-ExperimentNumber-1-MeteoCPC.csv
        Meteorological parameters and CPC measurements.
        From 60 s before expansion-start to 60 s after expansion-end
        Time (Sec): Seconds since expansion-start
        Pressure (hPa): Pressure in the cloud chamber
        Tgw_mean (degC): Mean value of gas temperature in the cloud chamber
        Tgw_std (degC): Standard deviation of gas temperature in the cloud chamber
        Tww_mean (degC): Mean value of wall temperature in the cloud chamber
        Tww_std (degC): Standard deviation of wall temperature in the cloud chamber
        TDew (degC), TFrost (degC): Dew and frost points measured by MBW
        CPC_Inters (#/cm3): Number concentration of interstitial particles [NOT YET INTERSTITIAL, CAN BE USED AS TOTAL]
        CPC_TotBot (#/cm3): Number concentration of all particles sampled at the bottom
        CPC_TotBot (#/cm3): Number concentration of all particles sampled at the top [CPC NOT WELL WORKING NOW]

    CampaignName-ExperimentNumber-2-OPC-*.csv
        Size distributions (dN/dlogDp, #/cm3) measured by each of the two Welas-Digital Sensors, as well as the one merged from two sensors.
        Optical equivalent diameters in um are given assuming the particles consist of water (refractive index = 1.33)
        For Welas01, dlogDp = [ 0.00913 0.00895 0.00855 0.00795 0.00775 0.00784 0.00767 0.00753 0.00740 0.00727 0.00715 0.00640 0.00561 0.00544 0.00526 0.00510
                                0.00495 0.00480 0.00468 0.00458 0.00508 0.00564 0.00567 0.00570 0.00572 0.00584 0.00596 0.00598 0.00600 0.00601 0.00601 0.00598
                                0.00595 0.00594 0.00585 0.00578 0.00556 0.00534 0.00540 0.00548 0.00558 0.00569 0.00582 0.00596 0.00611 0.00650 0.00678 0.00679
                                0.00686 0.00691 0.00695 0.00698 0.00700 0.00701 0.00702 0.00709 0.00717 0.00719 0.00718 0.00704 0.00689 0.00684 0.00676 0.00667
                                0.00659 0.00654 0.00651 0.00650 0.00662 0.00675 0.00678 0.00682 0.00687 0.00696 0.00706 0.00717 0.00722 0.00732 0.00757 0.00803
                                0.00849 0.00872 0.00885 0.00896 0.00906 0.00916 0.00882 0.00843 0.00840 0.00857 0.00985 0.01160 0.01272 0.01357 0.01405 0.01503
                                0.01625 0.01646 0.01625 0.01505 0.01588 0.01745 0.01672 0.01621 0.01565 0.01664 0.01871 0.01875 0.01832 0.01736 0.01598 0.01514
                                0.01557 0.01631 0.01465 0.01244 0.01190 0.01169 0.01047 0.00995 0.01070 0.01045 0.00961 0.00953 0.00949 0.00959 0.01064 0.01129
                                0.01062 0.01108 0.01155 0.01187 0.01254 0.01375 0.01490 0.01443 0.01451 0.01478 0.01363 0.01405 0.01467 0.01502 0.01445 0.01281
                                0.01341 0.01354 0.01304 0.01357 0.01346 0.01296 0.01202 0.01146 0.01217 0.01179 0.01132 0.01125 0.01143 0.01073 0.00975 0.00993
                                0.00983 0.00925 0.00878 0.00909 0.00954 0.00882 0.00869 0.00917 0.00894 0.00880 0.00901 0.00940 0.00894 0.00857 0.00902 0.00886
                                0.00906 0.00947 0.00953 0.00957 0.00954 0.00953 0.00991 0.01038 0.01053 0.01056 0.01072 0.01148 0.01173 0.01076 0.01077 0.01103
                                0.01066 0.01080 0.01038 0.01012 0.00975 0.00937 0.00990 0.00995 0.00955 0.00921 0.00881 0.00895 0.00876 0.00850 0.00857 0.00845
                                0.00814 0.00795 0.00805 0.00814 0.00799 0.00806 0.00818 0.00817 0.00799 0.00805 0.00846 0.00832 0.00836 0.00897 0.00881 0.00855
                                0.00894 0.00935 0.00960 0.00975 0.01019 0.01012 0.00963 0.00991 0.01023 0.00999 0.00980 0.00944 0.00910 0.00892 0.00887 0.00904
                                0.00901 0.00852 0.00799 0.00782 0.00784 0.00768 0.00751 0.00759 0.00732 0.00697 0.00719 0.00799 0.00697 0.00660 0.00781 0.00775]
        For Welas02, dlogDp = [ 0.00925 0.00907 0.00868 0.00803 0.00775 0.00781 0.00763 0.00749 0.00738 0.00728 0.00721 0.00718 0.00721 0.00726 0.00725 0.00722
                                0.00718 0.00700 0.00681 0.00676 0.00672 0.00669 0.00669 0.00672 0.00675 0.00679 0.00695 0.00713 0.00720 0.00727 0.00736 0.00742
                                0.00740 0.00744 0.00766 0.00798 0.00852 0.00902 0.00924 0.00941 0.00957 0.00974 0.00946 0.00912 0.00925 0.00893 0.01010 0.01223
                                0.01311 0.01380 0.01418 0.01506 0.01621 0.01635 0.01638 0.01524 0.01479 0.01622 0.01642 0.01594 0.01534 0.01611 0.01753 0.01793
                                0.01776 0.01673 0.01579 0.01516 0.01426 0.01474 0.01458 0.01289 0.01192 0.01161 0.01114 0.00993 0.00971 0.01058 0.01019 0.00921
                                0.00923 0.00925 0.00902 0.01001 0.01105 0.01046 0.00980 0.01042 0.01145 0.01196 0.01256 0.01384 0.01392 0.01345 0.01340 0.01339
                                0.01416 0.01419 0.01347 0.01283 0.01187 0.01133 0.01156 0.01267 0.01318 0.01232 0.01177 0.01149 0.01192 0.01212 0.01186 0.01186
                                0.01118 0.01082 0.01045 0.01067 0.01094 0.00997 0.00997 0.00999 0.00996 0.00976 0.00968 0.00982 0.00940 0.00994 0.01048 0.01020
                                0.01037 0.01029 0.01056 0.01067 0.01098 0.01147 0.01117 0.01149 0.01182 0.01184 0.01222 0.01225 0.01217 0.01176 0.01155 0.01197
                                0.01192 0.01180 0.01140 0.01104 0.01069 0.01034 0.01040 0.01027 0.00984 0.00929 0.00889 0.00898 0.00880 0.00868 0.00861 0.00841
                                0.00817 0.00795 0.00798 0.00795 0.00777 0.00761 0.00775 0.00779 0.00770 0.00776 0.00767 0.00773 0.00754 0.00748 0.00797 0.00789
                                0.00795 0.00821 0.00818 0.00786 0.00758 0.00834 0.00907 0.00891 0.00899 0.00917 0.00922 0.00931 0.00983 0.01026 0.00982 0.00960
                                0.00965 0.00946 0.00927 0.00925 0.00921 0.00918 0.00930 0.00903 0.00829 0.00813 0.00822 0.00795 0.00781 0.00785 0.00779 0.00712
                                0.00721 0.00775 0.00725 0.00713 0.00750 0.00765 0.00742 0.00706 0.00752 0.00768 0.00758 0.00811 0.00833 0.00808 0.00787 0.00819
                                0.00858 0.00855 0.00876 0.00918 0.00916 0.00900 0.00909 0.00923 0.00941 0.00944 0.00947 0.00935 0.00888 0.00865 0.00849 0.00853
                                0.00834 0.00781 0.00777 0.00765 0.00745 0.00732 0.00735 0.00749 0.00724 0.00712 0.00697 0.00734 0.00681 0.00604 0.00521 0.00400]
        For Merged, dlogDp = 0.00984

    CampaignName-ExperimentNumber-3-InitialPNSD-*.csv
        PNSD for SMPS, APS and the merged ones
        Size distributions (dN/dlogDp, #/cm3)     of dry particles at expansion-start, reported at the pressure and temperature in side the chamber at expansion-start.
        Size distributions (dN/dlogDp, #/std cm3) of dry particles at expansion-start, reported at the standard pressure and temperature (101325 Pa and 273.15 K).
        Volume equivalent diameters in um are given for merged PNSD assuming 
            ---- spherical particles with density of 1.77 g/cm3 for non-dust aerosol systems
            ---- dynamic shape factor of 1.50 and density of 2.60 g/cm3 for dust-containing aerosol systems
            ---- Those assumption is not good for some cases, and could be revised later with more detailed consideration.
        The number and range of sizes may change according to the instrument and aerosol properties.
