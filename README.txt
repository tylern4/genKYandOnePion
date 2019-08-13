

1. Compile: ./COMPILE_eg_ky

2. Adjust eg_config_test.txt file. 
 Its format is as follows.
 - first line must contains the string "KLambda" or "KSigma" or "Pi0P" or "PiN"
 - Next line contains Ebeam that should be >0 and less than 12GeV (Energy of the beam)
 - Next line contains two numbers Q2min and Q2max.
   This is the range in Q2, where the events will be generated.
 - Next line contains two numbers Wmin and Wmax.
   This is the range in W, where the events will be generated.
   Wmin more then zero (GeV); (if Wmin less than treshold, treshold value will be used)
   Wmax less then 4 GeV;
 - Next line contains the number of events to be generated.
 - Next line contains the name of the output file.
 - Next line contains the path name to the "data" directory.
   
3. Run it:
>  ./eg_ky eg_config_test.txt
  Parameter eg_config_test.txt is the name of the config. file

