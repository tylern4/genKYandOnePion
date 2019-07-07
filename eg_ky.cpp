//SYS LIBRARIES
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

//ROOT LIBRARIES
#include <TLorentzVector.h>

//#include "sigmaKY.h"
#include "sigmaValera.h"
#include "constants.h" 
#include "kinematics.h"
#include "evGenerator.h"


using namespace std;


int main(int argc, char **argv)
{


        int channel;
	string channelName, outputFileName;
	double Ebeam, Q2min, Q2max, Wmin, Wmax;
	int nEventMax;
        double jr, mr, gr, a12, a32, s12, onlyres;

        

        // parameters
	//channelName = "Pi0P"; //"KLambda" or "KSigma" or "Pi0P" or "PiN"
	channelName = "PiN"; //"KLambda" or "KSigma" or "Pi0P" or "PiN"
	Ebeam = 6.4;//>0
	//Ebeam = 10.6;//>0
	Q2min = 0.5;//more then 0.001
	Q2max = 7;//up to 12 GeV2;
	Wmin = 1.;// >=0
	Wmax = 2;//up to 4 GeV;
        nEventMax = 100000;
	//outputFileName = "lund_PiN_E10_f5t9_n2.lund";
	outputFileName = "lund_PiN_E6_f05t7_100_n4.lund";

	bool check_in_data=check_input_data(channelName,Ebeam,Q2min,Q2max,Wmin,Wmax,nEventMax);
	if (check_in_data==0) {return 0;}

	// initilize event generator
	evGenerator eg(channelName, Ebeam,  Q2min, Q2max, Wmin, Wmax);
	channel=num_chanel(channelName);

        cout << endl
	     << " eg_ky initialized for channel " << channelName << endl  
	     << " Beam energy: " << Ebeam << endl
	     << " Q2 range: " << Q2min <<"  "<< Q2max << ";  "  
	     << " W range: "  << Wmin  <<"  "<< Wmax << ";  " << endl
	     << "   Number of events " << nEventMax << endl; 
	     
 
        // 4-momenta of final electron, Kaon and Lambda, proton, pi- in LAB frame
	TLorentzVector Peini( TVector3(0.,0.,Ebeam), sqrt(Ebeam*Ebeam+0.005*0.005) );
	TLorentzVector Ppini( TVector3(0.,0.,0.),    0.9383);
	TLorentzVector Pefin, PK, PL, Ppfin, Ppim, Pgam;
	double Q2,W;
	// output  
	ofstream output(outputFileName.c_str()); 
	// Loop through events
        for (int i=0; i<nEventMax; i++) {	
		
	  // get event. 4-moeg_ky/menta of final state particle. 
	  // Values of Q2 and W are also returned.  
          eg.getEvent(Q2, W, Pefin, PK, PL);//, Ppfin, Ppim, Pgam);


          // output in lund format
	  // header
	  //output << "4 1 1 0 0 0 0 "
	output << "3 1 1 0 0 0 0 "
	         <<" "<< W <<" "<< Q2 <<" "<< getomega(Q2, W) 
		 << endl;
	  // electron
	  output 
	    << "1 -1 1 " << lundIdElectron << " 0 0 "
	    << Pefin.Px() <<" "<< Pefin.Py() <<" "<< Pefin.Pz() 
	    <<" "<< Pefin.E() <<" 0.0005"
	    <<" 0 0 0 "
	    << endl;

	if (channel==1){
	  // Kaon+
	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;

	  output 
	   // << "3 0 1 " << lundIdSigmaZero << " 0 0 "
	<< "3 0 1 " << lundIdLambda << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.115"
	    <<" 0 0 0 "
	    << endl;
	}

	if (channel==2){

	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;

	  output 
	    << "3 0 1 " << lundIdSigmaZero << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.192"
	    <<" 0 0 0 "
	    << endl;
	}

	if (channel==3){

	  output 
	    << "2 0 1 " << lundIdPiZero << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.134"
	    <<" 0 0 0 "
	    << endl;

	  output 
	    << "3 1 1 " << lundIdProton << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 0.9382"
	    <<" 0 0 0 "
	    << endl;
	}

	if (channel==4){

	  output 
	    << "2 1 1 " << lundIdPiPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.1395"
	    <<" 0 0 0 "
	    << endl;

	  output 
	    << "3 0 1 " << lundIdNeutron << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 0.939"
	    <<" 0 0 0 "
	    << endl;
	}



/*
	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;
	  if(channel==1 || channel==2 || channel==4){
	
	  // Lambda
	      //output 
	      //  << "3 0 1 " << lundIdLambda << " 0 0 "
	      //  << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	      //  <<" "<< PL.E() <<" 1.115"
	      //  <<" 0 0 0 "
	      //  << endl;
	  output 
	    << "3 1 1 " << lundIdProton << " 0 0 "
	    << Ppfin.Px() <<" "<< Ppfin.Py() <<" "<< Ppfin.Pz() 
	    <<" "<< Ppfin.E() <<" 0.9383"
	    <<" 0 0 0 "
	    << endl;

	
	  output 
	    << "4 -1 1 " << lundIdPiMinus << " 0 0 "
	    << Ppim.Px() <<" "<< Ppim.Py() <<" "<< Ppim.Pz() 
	    <<" "<< Ppim.E() <<" 0.1396"
	    <<" 0 0 0 "
	    << endl;
	


	  } else if(channel==3) {
	  // Sigma zero
	  output 
	    << "2 0 1 " << lundIdSigmaZero << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.192"
	    <<" 0 0 0 "
	    << endl;	  
	  }
	*/  
	  // test 
	  //TLorentzVector PfinSum = Pefin + PK + Ppfin + Ppim;
	  //cout << " Pfin " << PfinSum.E() <<" "<< PfinSum.Px() <<" "<< PfinSum.Py() <<" "<< PfinSum.Pz()
	  //     << endl; 
	  //cout << " miss 0 "  << (Peini + Ppini - Pefin - PK - Ppfin - Ppim).M() << endl;
	  //cout << " miss e "  << (Peini + Ppini - PK - Ppfin - Ppim).M() << endl;
	  //cout << " miss p "  << (Peini + Ppini - Pefin - PK - Ppim).M() << endl;
	  //cout << " miss K "  << (Peini + Ppini - Pefin - Ppfin - Ppim).M() << endl;
	  //cout << " miss pi " << (Peini + Ppini - Pefin - PK - Ppfin).M() << endl;
	  //cout << " Mass L "  << (Ppfin + Ppim).M() << endl;
	  
	  
	  if( int((i+1)/1000)*1000 == (i+1) | (i+1)==1) {
	    cout << " Event # " << i+1 << endl; 
	  }
	}
	output.close();
	

	return 0;
	
	
}
