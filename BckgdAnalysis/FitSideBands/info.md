This file contains how to use all codes in this location.

1. First run "AnalyzeMeanSigma.C". We will hardcode the range for different pT ranges in this code. After we pass hardcode pT ranges here, we do not have to do anything in subsequent codes. ZZZZ.csv file produced will take care of that.

This .C code will produce some png files which we can see to adjust .csv file so that we correctly get the range for extrapolation. So, output of this code/analysis is much required. Manual adjustment is always recommended at this stage but we can take care of that in subsequent steps too.

2. Second step is to run "ExpPol2Fit.C" that will do extrapolaton. This is the main area where you might want to adjust the range in .csv file and make sure the extrapolation is much better. This also produce the output image (not .csv file) and output image is useful to check the range for current fit function. The adjustment in .csv file is in the inva mass axis (ExcludeLeft , ExcludeRight for extrapolation).

After we think things are perfect, then we will deploy this for final analysis.

3. As the third step, we run "SignalBckgdHist.C" that will use .csv file and it produce output .root file that has signal and extrapolated background. This also subtract the bg from fg and extract signal (for eta peak) and fit using gauss function and integrate from 3 sigma left to 3 sigma to right to get yeild. These results are  kept in output .csv files. 

 
