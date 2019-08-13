# genKYandOnePion
EG is based on:
1) at Q2 < 5GeV2 interpolation of the available CLAS data 
2) at Q2 > 5GeV2 extrapolation of the available CLAS data that is based on Operator Product Expansion 

Requirments: Cern root.

Using:
1) install root (https://root.cern.ch/building-root) or type command: use root/6.10.02
2) git clone of the EG: git clone https://github.com/ValeriiKlimenko/genKYandOnePion/ 
3) Then type command: chmod +x COMPILE_eg_ky
4) To compile: ./COMPILE_eg_ky
5) to run exe file: ./eg_ky eg_config_test.txt

To set the parameters use eeg_config_test.txt
Adjust eg_config_test.txt file. 
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
	
The output will be a lund file with the name that you set in outputFileName variable.

Contact: valerii@jlab.org
