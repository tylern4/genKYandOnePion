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
5) to run exe file: ./eg_ky

To set the parameters use eg_ky.cpp file. There are parameters in eg_ky.cpp:

	channelName: it can be "KLambda" or "KSigma" or "Pi0P" or "PiN"
	Ebeam >0 and less than 12GeV (Energy of the beam)
	Q2min from 0.001 GeV2;  
	Q2max less then 11.99 GeV2;
	Wmin more then zero (GeV); (if Wmin less than treshold, treshold value will be used)
	Wmax less then 4 GeV;
        nEventMax >0;
	outputFileName should be something that ends on ".lund". An example: "pin_1.lund" 
	
The output will be a lund file with the name that you set in outputFileName variable.

Contact: valerii@jlab.org
