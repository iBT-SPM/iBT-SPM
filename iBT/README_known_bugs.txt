Known issues & bugs (last updated 2025-12-27)
----------------

* LI index value at max n for normal controls is not calculated or
displayed on summary plot if the deactivation
contrast is not also analysed. This is an unnecessary limitation
of the present code in iBrainTools_laterality_graph.m

* The first value of arrays of derived motion regressors or rejection
metrics that rely on knowing the position of the immediately
preceding scan (such as delayed regressors, DeltaOriginDistance
and FramewiseDisplacement) are currently set to zero by convention.
However, if the analysis "pick" range excludes some early scans after
preprocessing (as is typical), then we have available the
previous values and could therefore determine the true first value.
i.e. the metrics could be derived BEFORE the pick range
is applied, and then the pick selection applied to both the original
and derived metrics. Requires rearranging code in iBrainTools_stats.m
(between lines 860-977 of the version from 2020-05-24)

  
Other limitations
-----------------

* iBT has not been tested with SPM25 or later.

* The laterality section of the toolbox uses several functions 
requiring the Matlab Statistics toolbox (e.g. ttest, ttest2, jbtest).

* iBT_extract_TR.m requires MATLAB R2016b or later.

* As is the case for SPM generally, generation of FWEc and FDRc cluster-size
thresholded images is presently limited to one-tailed t-tests, and only
at 0.05 corrected.

* There are presently no iBT routines that utilise SPM12's default
spatial normalisation procedures. Our lab has moved to using
fMRIprep for preprocessing instead (see runinfo_fMRIPrep_norm.m).
iBT also continues to support the legacy routines that were the default
in SPM8. If you wish to perform SPM12's new default spatial normalisation
approach, please pre-process using SPM12's own batch engine, and then
utilize iBT for the statistical analysis and display.
