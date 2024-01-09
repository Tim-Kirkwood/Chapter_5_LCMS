# Chapter_5_LCMS
## Generate peakml files with retention time-aligned mass spectrometry peaks.
Raw LCMS data was copied from the LCMS machine onto a memory stick, and converted in MZML format with MSConvert from the proteowizard software suite.  This generated 65 MZML files, with file names follwoing the format ```host_bgc_replicate_runcode```.  Files were collected into one directory per strain (with each directory hiolding replicate samples for that strain) - directories were
```m1146_25, m1152_10, m1152_21, m145_10, m145_15, m145_21, m145_25, tk24_10, tk24_15, tk24_21, tk24_25```.  RT correction and PeakML file creation were then achieved via **pre_process.R**, which added directories containing individual PeakMl files, RT-corrected PeakMl files (for each replicate sample), and combined RT-coorected PeakMl files (for each strain) to each of the strain directories described above. 

## Annotation peak picking
This folder describes peak annotation with IPApy2 and subsequent search workflow for BGC-specific, confidently annotated peaks 

## Untargetted peak picking
This folder describes the search workflow for BGC-specific, confidently annotated peaks 
