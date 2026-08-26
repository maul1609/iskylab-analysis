"""Legacy entry point for the iSKYLAB batch analysis.

The original script duplicated the namelist-generation logic and depended on
obsolete BMM controls such as ``chamber_override``.  The maintained workflow is
now :mod:`dataAnalysis_new`; this wrapper is retained so old commands continue
to work without carrying a second stale implementation.
"""

from dataAnalysis_new import main


if __name__ == "__main__":
    main()
