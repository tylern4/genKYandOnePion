# genKYandOnePion
EG is based on:
1) at Q2 < 5GeV2 interpolation of the available CLAS data 
2) at Q2 > 5GeV2 extrapolation of the available CLAS data that is based on Operator Product Expansion 

Compile: ./COMPILE_eg_ky
run: ./eg_ky

To set the parameters use eg_ky.cpp file. There are parameters in eg_ky.cpp:

	channelName: it can be "KLambda" or "KSigma" or "Pi0P" or "PiN"
	Ebeam >0 and less than 12GeV (Energy of the beam)
	Q2min from 0.001 GeV2;  
	Q2max less then 11.99 GeV2;
	Wmin more then zero (GeV); (if Wmin less than treshold, treshold value will be used)
	Wmax less then 4 GeV;
        nEventMax >0;
