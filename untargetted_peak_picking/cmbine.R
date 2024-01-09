files_in = 'D:/all_mzmatch_data/negative/m145_samples/combined_RTcorr/21.peakml,D:/all_mzmatch_data/negative/m145_samples/combined_RTcorr/10.peakml,D:/all_mzmatch_data/negative/m145_samples/combined_RTcorr/15.peakml,D:/all_mzmatch_data/negative/m145_samples/combined_RTcorr/25.peakml'
file_out = 'D:/all_mzmatch_data/negative/m145_samples/combined_RTcorr/combined.peakml'
txt_file_out = 'D:/all_mzmatch_data/negative/m145_samples/combined_RTcorr/combined.txt'
mzmatch.ipeak.Combine(i = files_in, 
                      o = file_out,
                      rtwindow=60,
                      combination="set", 
                      ppm=50, 
                      v=T,#
                      nSlaves=2)

mzmatch.ipeak.sort.RelatedPeaks(i=file_out, v=T,
                                o="final_combined_related.peakml", 
                                basepeaks="final_combined_basepeaks.peakml", ppm=50, rtwindow=60)

mzmatch.ipeak.convert.ConvertToText (i=file_out, 
                                     o = txt_file_out)