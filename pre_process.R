install.packages ("remotes")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("https://github.com/andzajan/mzmatch.R.git",
                        build_opts=c("--no-multiarch"), INSTALL_opts=c("--no-test-load"))
library(mzmatch.R)
mzmatch.init(version.1=FALSE)

dir_base = 'D:/LCMS data processing'
samples = c('m1146_25', 
            'm1152_10',
            'm1152_21',
            'm145_10',
            'm145_15',
            'm145_21',
            'm145_25',
            'tk24_10',
            'tk24_15',
            'tk24_21',
            'tk24_25'
            )
#sample = "m145_10"
for (sample in samples) {
  setwd(paste(dir_base, '/', sample, sep = ''))

  tsv_file = paste(sample,"_sample_setup.tsv", sep = '')


  mzmatch.R.Setup(projectFolder= getwd(),
                  samplelist=tsv_file)
  
  xseto <- xcmsSet(sampleList$filenames, 
                   method='centWave', 
                   ppm=50,
                   peakwidth=c(10,60), 
                   snthresh=6,
                   mzdiff=0.01, 
                   integrate=1, 
                   prefilter=c(3,100),
                   
                   #from xcms params
                   noise = 0,
                   verbose.columns=TRUE,
                   fitgauss=FALSE, 
                   nSlaves=2)
  
  
  PeakML.xcms.write.SingleMeasurement(xset=xseto,
                                       filename=sampleList$outputfilenames, 
                                       ionisation="negative",
                                       ppm=50,
                                       addscans=0
                                       )
  
  
  mzmatch.R.Setup(projectFolder=getwd(), 
                  samplelist=tsv_file, 
                  outputfolder="peakml_RTcorr")
  xset2<-retcor(xseto, 
                method="obiwarp", 
                profStep=0.01, 
                center=3)
  PeakML.xcms.write.SingleMeasurement(xset=xset2, 
                                      filename=sampleList$outputfilenames, 
                                      ionisation="negative", 
                                      ppm=50,
                                      addscans=0, 
                                      )
  mzmatch.ipeak.Combine(sampleList=sampleList, 
                        rtwindow=60,
                        combination="set", 
                        ppm=50, 
                        v=T,#
                        nSlaves=2,
                        outputfolder="combined_RTcorr")
  inputfile = paste(getwd(), '/combined_RTcorr/', sample, '.peakml', sep='')
  output_file = paste(getwd(), '/combined_RTcorr/', sample, '.txt', sep='')
  mzmatch.ipeak.convert.ConvertToText (i = inputfile,
                                       o = output_file)
}