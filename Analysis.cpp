#define Analysis_cxx
// The class definition in Analysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Analysis.C")
// root> T->Process("Analysis.C","some options")
// root> T->Process("Analysis.C+")
//


#include "Analysis.h"

void Analysis::Begin(TTree * /*tree*/) {
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    /////////////////////Data//////////////////////

    Target = 1;
    BeData = 0;

    //////////////Target thickness////////////////
    TT10 = 0;
    TT20 = 0;
    TT30 = 0;
    TT38_6 = 0;
    TT43 = 0;//Target thickness 43 um
    TT45 = 1;//Target thickness 45 um
    TT47 = 0;//Target thickness 47 um
    TT49_99 = 0;//Target thickness 49.99 um
    TT56 = 1;//Target thickness 56 um for 12C from Geant4
    TT60 = 0;//Target thickness 60 um
    if(TT10){
        TThickness = 10;
    } else if(TT20) {
        TThickness = 20;
    } else if(TT30) {
        TThickness = 30;
    } else if(TT38_6) {
        TThickness = 38.6;
    } else if(TT43) {
        TThickness = 43;
    } else if(TT45) {
        TThickness = 45;
    } else if(TT47) {
        TThickness = 47;
    } else if(TT49_99) {
        TThickness = 49.99;
    } else if(TT56) {
        TThickness = 56.00;
    } else if(TT60) {
        TThickness = 60.78992;
    }

    /////////////////Si shift/////////////////////
    SiShift = 0;
    if(SiShift) {
        shift=0.0885;
    } else {
        shift=0;
    }

    ////////////////Beam offset////////////////////
    BeamOffset = 0;
    if(BeamOffset) {
        x_target = -3;
        y_target = 4;
    } else {
        x_target = 0;
        y_target = 0;
    }

    ////////////////Si Threshold////////////////////

    SiThresholdYu = 0.2;

    ///////////////Randomization////////////////////
    random = 1;

    ///////////////IC Channel Limits///////////////

    if(BeData){
        ICLowerBound = 620.;
        ICUpperBound = 1100.;
    }else{
        ICLowerBound = 1500.;
        ICUpperBound = 2200.;
    }
    

    ///////////////Energy Loss////////////////////

    if(BeData){
        //////////////////////SRIM///////////////////////
        // Beam_Al = new EnergyLoss("Be_in_Al.dat");
        // Beam_SiO2 = new EnergyLoss("Be_in_SiO2.dat");
        // Beam_B = new EnergyLoss("Be_in_B.dat");
        // Beam_Si = new EnergyLoss("Be_in_Si.dat");
        // Beam_P = new EnergyLoss("Be_in_P.dat");
        // Beam_Ag = new EnergyLoss("Be_in_Ag.dat");
        // Beam_D2 = new EnergyLoss("Be_in_D2_density0_201.dat");

        //////////////////////GEANT4///////////////////////
        Beam_Al = new EnergyLoss();
        Beam_SiO2 = new EnergyLoss();
        Beam_B = new EnergyLoss();
        Beam_Si = new EnergyLoss();
        Beam_P = new EnergyLoss();
        Beam_Ag = new EnergyLoss();
        Beam_D2 = new EnergyLoss();
        Beam_Al->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_Al_Energyloss.dat");
        Beam_SiO2->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_SiO2_Energyloss.dat");
        Beam_B->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_B_Energyloss.dat");
        Beam_Si->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_Si_Energyloss.dat");
        Beam_P->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_P_Energyloss.dat");
        Beam_Ag->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_Ag_Energyloss.dat");
        Beam_D2->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12Be_D2_Energyloss.dat");
    }else{

        //////////////////////SRIM///////////////////////
        // Beam_Al = new EnergyLoss("C_in_Al.dat");
        // Beam_SiO2 = new EnergyLoss("C_in_SiO2.dat");
        // Beam_B = new EnergyLoss("C_in_B.dat");
        // Beam_Si = new EnergyLoss("C_in_Si.dat");
        // Beam_P = new EnergyLoss("C_in_P.dat");
        // Beam_Ag = new EnergyLoss("C_in_Ag.dat");
        // Beam_D2 = new EnergyLoss("C_in_D2_density0_201.dat");

        //////////////////////GEANT4///////////////////////
        Beam_Al = new EnergyLoss();
        Beam_SiO2 = new EnergyLoss();
        Beam_B = new EnergyLoss();
        Beam_Si = new EnergyLoss();
        Beam_P = new EnergyLoss();
        Beam_Ag = new EnergyLoss();
        Beam_D2 = new EnergyLoss();
        Beam_Al->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_Al_Energyloss.dat");
        Beam_SiO2->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_SiO2_Energyloss.dat");
        Beam_B->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_B_Energyloss.dat");
        Beam_Si->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_Si_Energyloss.dat");
        Beam_P->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_P_Energyloss.dat");
        Beam_Ag->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_Ag_Energyloss.dat");
        Beam_D2->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_12C_D2_Energyloss.dat");
    }

    protonELB = new EnergyLoss();
    protonELA = new EnergyLoss();
    protonELD = new EnergyLoss();

    protonELB->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_1H_B_Energyloss.dat");
    protonELA->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_1H_Al_Energyloss.dat");
    protonELD->ReadBasicdEdx("Geant4EnergyLossTables/Geant4_1H_D2_Energyloss.dat");

    // protonELB = new EnergyLoss("Proton_Boron.dat");
    // protonELA = new EnergyLoss("Proton_Aluminum.dat");
    // protonELD = new EnergyLoss("Proton_DeuteriumTarget.dat");

    // protonELB->AddBackHigh(30.);
    // protonELA->AddBackHigh(30.);
    // protonELD->AddBackHigh(30.);


    // C_Al = new EnergyLoss("C_in_Al.dat");
    // C_SiO2 = new EnergyLoss("C_in_SiO2.dat");
    // C_B = new EnergyLoss("C_in_B.dat");
    // C_Si = new EnergyLoss("C_in_Si.dat");
    // C_P = new EnergyLoss("C_in_P.dat");

    // C_Al->AddBackError(1.e-4);
    // C_Al->UseGL128();
    // C_SiO2->AddBackError(1.e-4);
    // C_SiO2->UseGL128();
    // C_B->AddBackError(1.e-4);
    // C_B->UseGL128();
    // C_Si->AddBackError(1.e-4);
    // C_Si->UseGL128();
    // C_P->AddBackError(1.e-4);
    // C_P->UseGL128();


    // C_Ag = new EnergyLoss("C_in_Ag.dat");
    // C_D2 = new EnergyLoss("C_in_D2_density0_201.dat");

    // C_Ag->AddBackError(1.e-4);
    // C_Ag->UseGL128();

    // Be_Al = new EnergyLoss("Be_in_Al.dat");
    // Be_SiO2 = new EnergyLoss("Be_in_SiO2.dat");
    // Be_B = new EnergyLoss("Be_in_B.dat");
    // Be_Si = new EnergyLoss("Be_in_Si.dat");
    // Be_P = new EnergyLoss("Be_in_P.dat");

    // Be_Al->AddBackError(1.e-4);
    // Be_Al->UseGL128();
    // Be_SiO2->AddBackError(1.e-4);
    // Be_SiO2->UseGL128();
    // Be_B->AddBackError(1.e-4);
    // Be_B->UseGL128();
    // Be_Si->AddBackError(1.e-4);
    // Be_Si->UseGL128();
    // Be_P->AddBackError(1.e-4);
    // Be_P->UseGL128();

    // Be_Ag = new EnergyLoss("Be_in_Ag.dat");
    // Be_D2 = new EnergyLoss("Be_in_D2_density0_201.dat");

    // Be_Ag->AddBackError(1.e-4);
    // Be_Ag->UseGL128();

    ////////Sd Deadlayer thicknesses//////////
    dAg = 4.6405042/1000.;

    //Deadlayer D1
    dSdAl1=1.5/1000.;
    dSdSiO2=3.5/1000.;
    dSdAl2=0.3/1000.;
    dSdB=0.5/1000.;

    //Active Si layer
    dSd1Si=61./1000.;
    dSd2Si=493./1000.;

    //Deadlayer D2
    dSdP=0.5/1000.;
    dSdAl3=0.3/1000.;

    //read in the geometry file to calculate the angles
    std::ifstream geometry;
    geometry.open ("geometry_s1506.txt"); //open the geometry file; change the path and name accordingly
    if(!geometry.is_open()){
        std::cout << " No Geometry file found " << std::endl;
        return;
    }

    std::string read_geometry;
    std::istringstream iss;
    std::string name, dummy;

    while(getline(geometry, read_geometry)) { //in the calibration file named geometry start reading the lines
        if(read_geometry.find("YD_DISTANCE", 0) != std::string::npos) { //if you find the "YD_DISTANCE" before the end of the file
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Ydz;
        }

        if(read_geometry.find("YD_INNER_RADIUS", 0) != std::string::npos) {
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Ydr0;
        }

        if(read_geometry.find("YD_OUTER_RADIUS", 0) != std::string::npos) {
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Ydr1;
        }

        if(read_geometry.find("YU_DISTANCE", 0) != std::string::npos) {
            iss.clear();
            iss.str ( read_geometry );
            iss >> name >> dummy >> Yuz;
        }

        if(read_geometry.find("SD1_DISTANCE", 0) != std::string::npos) {
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Sd1z;
        }

        if(read_geometry.find("SD2_DISTANCE", 0) != std::string::npos) {
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Sd2z;
        }

        if(read_geometry.find("SD_INNER_RADIUS", 0) != std::string::npos) {
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Sdr0;
        }

        if(read_geometry.find("SD_OUTER_RADIUS", 0) != std::string::npos) {
            iss.clear();
            iss.str(read_geometry);
            iss >> name >> dummy >> Sdr1;
        }
    }
    //end of the geometry file

    ////////Output file//////////

    if(BeData){
        f_cut_name = "../BeCutIC2.root";
        f_cut = TFile::Open(f_cut_name); 
        pidcut = (TCutG*) f_cut->Get("BeCut_NewCalSd1rSd2r");//PID cut for Be
        if(Target){
            if(random) {
                // f_out_name = Form(path+"/Analysis/BeamOffset/beryllium/Be_pedestal_TRIUMF_DL_BeamOffset_%.0f_%.0f_Random_Shift_%.4f_CutIC_TargetDistance%.2fmm_TargetThickness%.2fum_Yu_Mod_UnevenBinning_Sd1rAlphaCal_Sd2rNewInBeamCal_NewSectorGeometry_ElasticScattering.root",x_target,y_target,shift,0-Yuz,TThickness);
                f_out_name = "Be_output_random_T_Yd.root";
            }else{
                // f_out_name = Form(path+"/Analysis/BeamOffset/beryllium/Be_pedestal_TRIUMF_DL_BeamOffset_%.0f_%.0f_Shift_%.4f_CutIC_TargetDistance%.2fmm_TargetThickness%.2fum_Yu_Mod_UnevenBinning_Sd1rAlphaCal_Sd2rNewInBeamCal_NewSectorGeometry_ElasticScattering.root",x_target,y_target,shift,0-Yuz,TThickness);
                f_out_name = "Be_output_T_AlphaCheck.root";
            }
        }else{
            if(random) {
                f_out_name = "Be_output_random_NT_AlphaCheck.root";
            }else{
                f_out_name = "Be_output_NT_AlphaCheck.root";
            }
        }
    }else{
        f_cut_name = "../CCalCut.root";
        f_cut = TFile::Open(f_cut_name); 
        pidcut = (TCutG*) f_cut->Get("CcutCalSd1rSd2rFull");//PID cut for C
        if(Target){
            if(random) {
                f_out_name = "C_output_random_T.root";
            }else{
                f_out_name = "C_output_T.root";
            }
        }else{
            if(random) {
                f_out_name = "C_output_random_NT.root";
            }else{
                f_out_name = "C_output_NT.root";
            }
        }
    }

    f_out = new TFile(f_out_name, "RECREATE");
    // g_out = new TFile("C_YuSiEnergy.root", "RECREATE");

    f_cutProton = TFile::Open("../ProtonCut.root"); //PID cut for tritons (Yd)
	pidcutProton = (TCutG*) f_cutProton->Get("protonCut");

    f_cutTriton = TFile::Open("../TritonCut.root"); //PID cut for tritons (Yd)
	pidcutTriton = (TCutG*) f_cutTriton->Get("tritonCut");

    f_cutDeuteron = TFile::Open("../DeuteronCut.root"); //PID cut for Deuterons (Yd)
	pidcutDeuteron = (TCutG*) f_cutDeuteron->Get("deuteronCut");

    f_cutElastic = TFile::Open("../ElasticCut.root"); //PID cut for Elastics (Sd1r vs Sd2r)
	pidcutElastic = (TCutG*) f_cutElastic->Get("CSd1rCSd2rElastic");

    f_cutHorizontal1 = TFile::Open("../HorizontalCut1.root"); //PID cut for Elastics (Sd1r vs Sd2r)
	pidcutHorizontal1 = (TCutG*) f_cutHorizontal1->Get("horizontalCut1");

    f_cutHorizontal2 = TFile::Open("../HorizontalCut2.root"); //PID cut for Elastics (Sd1r vs Sd2r)
	pidcutHorizontal2 = (TCutG*) f_cutHorizontal2->Get("horizontalCut2");

    f_cutHorizontal3 = TFile::Open("../HorizontalCut3.root"); //PID cut for Elastics (Sd1r vs Sd2r)
	pidcutHorizontal3 = (TCutG*) f_cutHorizontal3->Get("horizontalCut3");

    //definition of histograms
    //calculate the variable bins for the Yu detector to plot the angles
    width = (Ydr1 - Ydr0)/16;

	for(int i=0; i<17;i++){
		Ydbins[i]=TMath::RadToDeg()*TMath::ATan((Ydr0+i*width)/Ydz);
	}
    
	// for(int k=16; k>=0; k--){
	// 	int a = TMath::Abs(k-16);
	// 	Ydbins[a]=Ydbintemp[k];
	// }


    for(int i = 0; i < 17; i++) {
        Yubins[i] = 180 - TMath::RadToDeg()*TMath::ATan((50 + (16 - i)*width)/(0 - Yuz));
        //std::cout << Yubins[i] << std::endl;
    }
    for(int i = 0; i < 17; i++) {
        YubinsRad[i] = TMath::DegToRad()*180 - TMath::ATan((50 + (16 - i)*width)/(0 - Yuz));
        //std::cout << Yubins[i] << std::endl;
    }

    widthSd = (Sdr1 - Sdr0)/24;
    for(int i = 0; i < 25; i++) {
        Sd1Bins[i] = TMath::RadToDeg()*TMath::ATan((Sdr0 + i*widthSd)/Sd1z);
        //std::cout << Sd1Bins[i] << std::endl;
    }
    for(int i = 0; i < 24; i++) {
        Sd1Angles[i] = TMath::RadToDeg()*TMath::ATan((Sdr0 + (i + 0.5)*widthSd)/Sd1z);
        //std::cout << Sd1Angles[i] << std::endl;
    }

    for(int i = 0; i < 25; i++) {
        Sd2Bins[i] = TMath::RadToDeg()*TMath::ATan((Sdr0 + i*widthSd)/Sd2z);
        //std::cout << Sd2Bins[i] << std::endl;
    }

    for(int i = 0; i < 24; i++) {
        Sd2Angles[i] = TMath::RadToDeg()*TMath::ATan((Sdr0 + (i + 0.5)*widthSd)/Sd2z);
        //std::cout << Sd2Angles[i] << std::endl;
    }

    //Solid angle calculation
    int n=1;
    int distLength=(16/n) + 1, omegaLength=(16/n);
    double dist[distLength], omega[omegaLength];

    for(int i = 0; i < distLength; i++) {
        dist[i] = Ydr0 + n*i*width;
    }
    for(int i = 0; i < omegaLength; i++){
        omega[i] = TMath::DegToRad()*42*(0 - Yuz)*((1/(pow(TMath::Sqrt(dist[i]), 2) + Yuz*Yuz)) - (1/(pow(TMath::Sqrt(dist[i + 1]), 2) + Yuz*Yuz)));
    }
    double solidAngleBins[distLength];
    for(int i = 0; i < distLength; i++){
        dist[i] = Ydr0 + n*i*width;
        solidAngleBins[i] = 180 - TMath::RadToDeg()*TMath::ATan((Ydr0 + (distLength - 1 - i)*n*width)/(0 - Yuz));
        //std::cout << solidAngleBins[i] << std::endl;
    }

    ////////////////
    // Histograms //
    ////////////////

    hSd1rSd2r = new TH2D("hSd1rSd2r","Raw Sd1r vs Sd2r",1024,0,4096,1024,0,4096);
    hCSd1rSd2r = new TH2D("hCSd1rSd2r","Calibrated Sd1r vs Sd2r",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r
    hCSd1rSd2rElastic = new TH2D("hCSd1rSd2rElastic","Calibrated Sd1r vs Sd2r elastic gate",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r elastic gate
    hCSd1rSd2rProton_Yd = new TH2D("hCSd1rSd2rProton_Yd","Calibrated Sd1r vs Sd2r gated by Yd ProtonCut",1500,0,150,600,0,60);
    hCSd1rSd2rDeuteron_Yd = new TH2D("hCSd1rSd2rDeuteron_Yd","Calibrated Sd1r vs Sd2r gated by Yd DeuteronCut",1500,0,150,600,0,60); // calibrated Sd1r vs Sd2r gated by Yd DeuteronCut
    hCSd1rSd2rTriton_Yd = new TH2D("hCSd1rSd2rTriton_Yd","Calibrated Sd1r vs Sd2r gated by Yd TritonCut",1500,0,150,600,0,60);

    hCSd1rSd2rHorizontal1 = new TH2D("hCSd1rSd2rHorizontal1","Calibrated Sd1r vs Sd2r gated by HorizontalCut1",1500,0,150,600,0,60);
    hCSd1rSd2rHorizontal2 = new TH2D("hCSd1rSd2rHorizontal2","Calibrated Sd1r vs Sd2r gated by HorizontalCut2",1500,0,150,600,0,60);
    hCSd1rSd2rHorizontal3 = new TH2D("hCSd1rSd2rHorizontal3","Calibrated Sd1r vs Sd2r gated by HorizontalCut3",1500,0,150,600,0,60);

    hCSd1rSd2rIC = new TH2D("hCSd1rSd2rIC","Cal Sd1r vs Sd2r, gated by IC",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r and gated by IC
    hCSd1rSd2rICCut = new TH2D("hCSd1rSd2rICCut","Cal Sd1r vs Sd2r, gated by IC gated by the channel corr",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r gated by IC and required to be
    hCSd1rSd2rYdIC = new TH2D("hCSd1rSd2rYdIC","Cal Sd1r vs Sd2r, gated by Yd & IC",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r and gated by Yd & IC
    hCSd1rSd2rYdICCut = new TH2D("hCSd1rSd2rYdICCut","Cal Sd1r vs Sd2r, gated by Yd & IC & chann corr",1500,0,150,600,0,40); // calibrated Sd1r vs Sd2r gated by Yd, IC and required to be correlated with each other

    hYuMul = new TH1D("hYuMul","hYuMul",10,0,10);

    hYuEn = new TH1D("hYuEn","hYuEn",2500,0,10); //energy singles in Yu
    hYuEnT = new TH1D("hYuEnT","hYuEn with target energy loss",2500,0,10);
    hYuEn1 = new TH1D("hYuEn1","hYuEn1",2500,0,10);

    hYuEnPID = new TH1D("hYuEnPID","hYuEnPID",2400,0,2.4);//6000,0,6; 2400,0,2.4

    hYuEnIC = new TH1D("hYuEnIC","hYuEnIC",2400,0,2.4);
    hYuEnIC1 = new TH1D("hYuEnIC1","hYuEnIC1",2400,0,2.4);

    hYuAn = new TH2D("hYuAn","YuE vs Angle",16,Yubins,250,0,10); //Energy vs angle in Yu no gates
    hYuAnT = new TH2D("hYuAnT","YuE vs Angle with target energy loss",16,Yubins,250,0,10);
    hYuAn1 = new TH2D("hYuAn1","YuE vs Angle",16,Yubins,250,0,10);

    hYuAnIC = new TH2D("hYuAnIC","YuE vs Angle with a gate on the IC",16,Yubins,6000,0,6); //Energy vs angle in Yu with a gate on IC//2400,0,2.4 for 12Be
    hYuAnICT = new TH2D("hYuAnICT","YuE vs Angle with a gate on the IC and with target energy loss",16,Yubins,240,0,2.4);//240,0,2.4
    hYuAnIC1 = new TH2D("hYuAnIC1","YuE vs Angle with a gate on the IC",16,Yubins,240,0,2.4);//240,0,2.4

    hYuAnPID_SiEn = new TH2D("hYuAnPID_SiEn","Yu Silicon Energy vs Angle with an IC and PID gate",16,Yubins,6000,0,6);

    hYuAnPID = new TH2D("hYuAnPID","YuE vs Angle with an IC and PID gate",16,Yubins,6000,0,6); //Energy vs angle in Yu with a gate on IC and PID//240,0,2.4
    hYuAnPIDRad = new TH2D("hYuAnPIDRad","YuE vs Angle (in rad) with an IC and PID gate",16,YubinsRad,240,0,2.4);
    hYuAnPIDT = new TH2D("hYuAnPIDT","YuE vs Angle with an IC and PID gate and with target energy loss",16,Yubins,240,0,2.4);//240,0,2.4
    hYuAnPID1 = new TH2D("hYuAnPID1","YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4);//240,0,2.4
    hYuAnPIDYd = new TH2D("hYuAnPIDYd","YuE vs Angle with IC, Yd, Su, and PID gates",16,Yubins,240,0,2.4); //Energy vs angle in Yu with a gate on IC and PID//240,0,2.4
    hYuAnPID0_90 = new TH2D("hYuAnPID0_90","YuE vs Angle with an IC and PID gate, 0-90 degrees",16,Yubins,240,0,2.4);
    hYuAnPID90_180 = new TH2D("hYuAnPID90_180","YuE vs Angle with an IC and PID gate, 90-180 degrees",16,Yubins,240,0,2.4);
    hYuAnPID180_270 = new TH2D("hYuAnPID180_270","YuE vs Angle with an IC and PID gate, 180-270 degrees",16,Yubins,240,0,2.4);
    hYuAnPID270_360 = new TH2D("hYuAnPID270_360","YuE vs Angle with an IC and PID gate, 270-360 degrees",16,Yubins,240,0,2.4);
    hYuAnPIDRs0_7 = new TH2D("hYuAnPIDRs0_7","YuE vs Angle with an IC and PID gate, rings 0 to 7",16,Yubins,240,0,2.4);
    hYuAnPIDRs8_15 = new TH2D("hYuAnPIDRs12_15","YuE vs Angle with an IC and PID gate, rings 8 to 15",16,Yubins,240,0,2.4);
    hSd1An = new TH2D("hSd1An","Sd1E vs Angle",24,Sd1Bins,250,0,25);
    hSd2An = new TH2D("hSd2An","Sd2E vs Angle",24,Sd2Bins,1200,0,120);

    hCSd1rAngleElastic = new TH1D("hCSd1rAngleElastic","Sd1 Angle Elastic Cut",40,0,4);
    hCSd2rAngleElastic = new TH1D("hCSd2rAngleElastic","Sd2 Angle Elastic Cut",40,0,4);

    hYdAngleElastic = new TH1D("hYdAngleElastic","Yd Angle Elastic Cut",450,0,90);

    hSdReconAn = new TH2D("hSdReconAn","Reconstructed Sd Energy vs Angle",24,Sd1Bins,1200,0,120);
    hSdReconT = new TH1D("hSdReconT","Reconstructed Sd Energy With Target",1200,0,120);

    hSdReconAfterSilverAn = new TH2D("hSdReconAfterSilverAn","Reconstructed Sd Energy after silver vs Angle",24,Sd1Bins,1200,0,120);
    hSdReconAfterSilver = new TH1D("hSdReconAfterSilver","Reconstructed Sd Energy after silver",1200,0,120);
    hSd1NoDL1An = new TH2D("hSd1NoDL1An","Reconstructed Sd1 Energy with no DL1 vs Angle",24,Sd1Bins,1200,0,120);
    hSd1NoDL1 = new TH1D("hSd1NoDL1","Reconstructed Sd1 Energy with no DL1",1200,0,120);
    hSd1AfterActiveLayerAn = new TH2D("hSd1AfterActiveLayerAn","Reconstructed Beam Energy after Sd1 Active Layer vs Angle",24,Sd1Bins,1200,0,120);
    hSd1AfterActiveLayer = new TH1D("hSd1AfterActiveLayer","Reconstructed Beam Energy after Sd1 Active Layer",1200,0,120);

    hSdReconSd1An = new TH2D("hSdReconSd1An","Reconstructed Sd Energy vs Angle from Sd1rEn",24,Sd1Bins,1200,0,120);
    hSdReconSd1 = new TH1D("hSdReconSd1","Reconstructed Sd Energy from Sd1rEn",1200,0,120);

    hYdCsI1 = new TH2D("hYdCsI1","Yd vs CsI2",1000,10,60,250,0,10); // Yd vs CsI to see the elastic scattering
    hYdCsI2 = new TH2D("hYdCsI2","Yd vs CsI2",1000,0,50,250,0,10);
    hYdCsI1Elastic = new TH2D("hYdCsI1Elastic","Yd vs CsI1 with Elastic Gate",1000,10,60,250,0,10);
    hYdAn = new TH2D("hYdAn","YdE vs Angle",16,Ydbins,250,0,10);
    hYdAnPIDProton = new TH2D("hYdAnPIDProton","YdE vs Angle with the Proton PID cut",16,Ydbins,250,0,10);
    hYdAnPIDDeuteron = new TH2D("hYdAnPIDDeuteron","YdE vs Angle with the Deuteron PID cut",16,Ydbins,250,0,10);
    hYdAnPIDTriton = new TH2D("hYdAnPIDTriton","YdE vs Angle with the Triton PID cut",16,Ydbins,250,0,10);
    hYdAnPIDElastic = new TH2D("hYdAnPIDElastic","YdE vs Angle with the Elastic Gate",16,Ydbins,250,0,10);
    for(int i = 0; i < 16; i++) {
        hYdElastic[i] = new TH1D(Form("hYdElasticR%d",i),Form("Yd Energy, ring %d",i),250,0,10);
    }

    hYdCsI1_Tot = new TH2D("hYdCsI1_Tot","Yd vs Yd + CsI1",1000,10,60,250,0,10);
    hCsIAnPIDDeuteron = new TH2D("hCsIAnPIDDeuteron","CsI1 vs Angle with the Deuteron PID cut",16,Ydbins,1200,0,60);
    hYdCsI1Tot_AnPIDDeuteron = new TH2D("hYdCsI1Tot_AnPIDDeuteron","YdE +CsI1 vs Angle with the Deuteron PID cut",16,Ydbins,1200,0,60);


    // Q values
    hEx = new TH1D("hEx","Excitation energy",500,-5,5);
    hqvalvsEx = new TH2D("hqvalvsEx","Q value vs. Excitation energy ",500,-5,5,500,-5,5);
    float qval_i = -5.;
    float qval_f = -2.;
    float qval_bins = 50.;
    if(BeData){
        hQval_cut = new TH1D("hQval_cut",Form("Q values with the cut and Si threshold %0.2f MeV",SiThresholdYu),qval_bins,qval_i,qval_f);//qval_bins,qval_i,qval_f for Be
        hQval = new TH1D("hQval","Q values",qval_bins,qval_i,qval_f); 
    }else{
        hQval_cut = new TH1D("hQval_cut",Form("Q values with the cut and Si threshold %0.2f MeV",SiThresholdYu),150,-7,4);//qval_bins,qval_i,qval_f for Be
        hQval = new TH1D("hQval","Q values",150,-7,4);//100,-5,0
    }
    
    // hQval->GetYaxis()->SetTitle(Form("Counts/%.0f keV",(qval_f-qval_i)/qval_bins));
    hQval0_90 = new TH1D("hQval0_90","Q values, 0-90 degrees",800,-5,0);
    hQval90_180 = new TH1D("hQval90_180","Q values, 90-180 degrees",800,-5,0);
    hQval180_270 = new TH1D("hQval180_270","Q values, 180-270 degrees",800,-5,0);
    hQval270_360 = new TH1D("hQval270_360","Q values, 270-360 degrees",800,-5,0);
    hQval0_7 = new TH1D("hQval0_7","Q values (Rings 0 to 7, angles 134.5 to 149.5)",800,-5,0);
    hQval8_15 = new TH1D("hQval8_15","Q values (Rings 8 to 15, angles 123.5 to 134.5)",800,-5,0);
    hQvalT = new TH1D("hQvalT","Q values with target energy loss",400,-5,0);
    hQval1 = new TH1D("hQval1","Q values",400,-5,0);
    hQvalAn = new TH2D("hQvalAn","Q values vs Angle with energy loss",16,Yubins,400,-5,0);
    hQvalAnT = new TH2D("hQvalAnT","Q values vs Angle with target energy loss",16,Yubins,400,-5,0);
    hQvalYd = new TH1D("hQvalYd","Q values gated by Yd and Su",400,-5,0);

    hYuEnM = new TH2D("hYuEnM","YuE1 vs YuE2 for multiplicity 2",2400,0,2.4,2400,0,2.4);
    hYuEnMICPID = new TH2D("hYuEnMICPID","YuE1 vs YuE2 for multiplicity 2 with IC and PID gates",2400,0,2.4,2400,0,2.4);//2400,0,2.4

    hCSd1rEn = new TH1D("hCSd1rEn","Calibrated SD1r Energy",500,0,50);
    for (int i=0;i<24;i++){
        hCSd1rEnRing[i] = new TH1D(Form("hCSd1rEnRing%d",i),Form("Sd1r Calbrated energy, ring %d",i),500,0,50);
    }
    hCSd1r_0_1_2_En = new TH1D("hCSd1r_0_1_2_En","Calibrated SD1r Energy (Rings 0, 1, and 2)",500,0,50);

    hCSd1sEn = new TH1D("hCSd1sEn","Calibrated SD1s Energy",500,0,50);

    hCSd2rEn = new TH1D("hCSd2rEn","Calibrated SD2r Energy",500,0,120);
    for (int i=0;i<24;i++){
        hCSd2rEnRing[i] = new TH1D(Form("hCSd2rEnRing%d",i),Form("Sd2r Calbrated energy, ring %d",i),500,0,120);
    }
    hCSd2r_0_1_2_En = new TH1D("hCSd2r_0_1_2_En","Calibrated SD2r Energy (Rings 0, 1, and 2)",500,0,120);

    hCSd2sEn = new TH1D("hCSd2sEn","Calibrated SD2s Energy",500,0,120);

    hCSdr = new TH1D("hCSdrEn","Calibrated Sd1r+Sd2r Energy",500,0,120);

    for(int i = 0; i < 16; i++) {
        hQvalR[i] = new TH1D(Form("hQvalR%d",i),Form("Q values, ring %d",i),800,-5,0);
    }

    for(int i = 0;i < 8; i++) {
        hQval2R[i] = new TH1D(Form("hQval2R_%d_%d",2*i,2*i+1),Form("Q values, rings %d and %d",2*i,2*i+1),800,-5,0);
    }

    for(int i = 0; i < 5; i++) {
        hQval3R[i] = new TH1D(Form("hQval3R_%d_%d_%d",3*i,3*i+1,3*i+2),Form("Q values, rings %d , %d and %d",3*i,3*i+1,3*i+2),800,-5,0);
    }

    for (int i=0;i<4;i++){
        hQval4R[i] = new TH1D(Form("hQval4R_%d_%d_%d_%d",4*i,4*i+1,4*i+2,4*i+3),Form("Q values, rings %d , %d , %d and %d",4*i,4*i+1,4*i+2,4*i+3),800,-5,0);
    }

    for (int i=0;i<8;i++){
        hQvalS[i] = new TH1D(Form("hQvalS%d",i),Form("Q values, sector %d",i),800,-5,0);
    }

    for(int i=0; i<8; i++){
        std::string namehYuAnPIDSec = Form("hYuAnPIDSec_%i",i);
        hYuAnPIDSec[i] = new TH2D(namehYuAnPIDSec.c_str(),"YuE vs Angle with an IC and PID gate",16,Yubins,240,0,2.4);//240,0,2.4
    }

    hAngle_1 = new TH1D("hAngle_1", "Counts vs. Angle for the State 1; Angle [deg]; Counts",omegaLength,solidAngleBins);
    hAngle_1->Sumw2();

    hAngle_2 = new TH1D("hAngle_2", "Counts vs. Angle for the State 2; Angle [deg]; Counts",omegaLength,solidAngleBins);
    hAngle_2->Sumw2();

    hAngle_3 = new TH1D("hAngle_3", "Counts vs. Angle for the State 3; Angle [deg]; Counts",omegaLength,solidAngleBins);
    hAngle_3->Sumw2();

    hAngle_4 = new TH1D("hAngle_4", "Counts vs. Angle for the State 4; Angle [deg]; Counts",omegaLength,solidAngleBins);
    hAngle_4->Sumw2();

    ev_num = 0;

    //variable calculation for the YY1 detectors used in angle calculations
    YChWidth = (Ydr1 - Ydr0)/16.;
    YChMiddle = YChWidth/2;
    yuM = 0; ydM = 0;
    beam.SetX(0); beam.SetY(0); beam.SetZ(1);

    //variable calculation for the S3 detectors used in angle calculations
    SdChWidth = (Sdr1 - Sdr0)/24.;
    SdChMiddle = SdChWidth/2;
    SdM = 0;

    //variable definition for the Q value calculations
    amu = 931.5; // atomic mass unit in MeV
    massEjec = 938.28; //mass of the proton in MeV/c2
    

    if(BeData){
        //Be
        if(Target){
            if(TT43) {
                kBeam = 112.37;
            } else if(TT45) {
                kBeam = 112.35;
            } else if(TT47) {
                kBeam = 112.33;
            } else if(TT49_99) {
                kBeam = 112.31;
            } else if(TT60) {
                kBeam = 112.21;
            } else if(TT10) {
                kBeam = 112.66;
            }
        }else{
            kBeam = 112.75;
        }
    }else{
        //C
        if(Target){
            if(TT43) {
                kBeam = 110.36;
            } else if(TT49_99) {
                kBeam = 110.22;
            }else if(TT56){//Geant4 12C
                kBeam = 110.10;
            }else if(TT30){
                kBeam = 110.62;
            }else if(TT20){
                kBeam = 110.82;
            }else if(TT10){
                kBeam = 111.02;
            }else if(TT38_6){
                kBeam = 110.45;
            }
        }else{
            kBeam = 111.22;
        }
    }

    mAg = 107.87*amu;//mass of the silver foil

    mbeam = 12 * amu;  //mass of the beam (12Be or 12C) in MeV
    mrecoil = 13 * amu;  //mass of the recoil (13Be or 13C) in MeV
    mejec = 1 * amu; //mass of the proton

    //////////////////Calculations for finding the energy of the beam scattered in the silver////////////////////
    AScat = pow(mbeam/(mbeam + mAg),2);
    BScat = pow(mAg/mbeam,2);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << TThickness << " " << kBeam << " " << 0 - Yuz << std::endl;

    phi[0]=90.0;
    for(int i = 1; i < 8; i++) {
        phi[i] = phi[0] - 45*i;
        //std::cout << phi[i] << std::endl;
    }

    for(int i = 0; i < 16; i++) {
        r[i] = 50 + ((129. - 50.)/16)*(0.5 + i);
        //std::cout << r[i] << std::endl;
    }

    //Beam offset
    for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 16; j++) {
            double x_yu = r[j]*cos(phi[i]*M_PI/180.);
            double y_yu = r[j]*sin(phi[i]*M_PI/180.);

            TVector3 vec1(0, 0, 0 - Yuz);
            TVector3 vec2(x_yu - x_target, y_yu - y_target, Yuz);
            double angle = vec1.Angle(vec2);
            // std::cout << i << '\t' << j << '\t' << angle*180./M_PI << std::endl;

            DetAngle[i][j] = angle*180./M_PI;
        }
    }

    seed = 123;
    ran = new TRandom3(seed);

    Q13Be0 = -2.649;

    // Read in new calibration for Sd1 and Sd2

    int ring, channel_csi1, sector_csi1;
    double gain, offset;

    if(BeData){
        std::ifstream sd1Calibration("Sd1rAlphaCalibration.txt");
        // std::ifstream sd1Calibration("Sd1r_5point_Calibration.txt");
        assert(sd1Calibration.is_open());
        while(sd1Calibration >> ring >> gain >> offset) {
        Sd1rRingCalibration[ring] = std::make_pair(gain, offset);
        }
        sd1Calibration.close();
    }else{
        std::ifstream sd1Calibration("Sd1rInBeamCalibration.txt");
        assert(sd1Calibration.is_open());
        while(sd1Calibration >> ring >> gain >> offset) {
        Sd1rRingCalibration[ring] = std::make_pair(gain, offset);
        }
        sd1Calibration.close();
    }
    

    std::ifstream csi1Calibration("calibration_CsI1_500x_jerome.txt");
    assert(csi1Calibration.is_open());
    while(csi1Calibration >> channel_csi1 >> ring >> offset >> gain) {
        sector_csi1 = channel_csi1-32;   //csi1 channel number starts from 32 in the calibration file
        // Csi1RingSectorCalibration[16*sector_csi1+ring] = std::make_pair(gain, offset);
        Csi1RingSectorCalibration[sector_csi1][ring] = std::make_pair(gain, offset);
        // std::cout << channel_csi1 << "\t" << sector_csi1 << "\t" << ring << "\t" << Csi1RingSectorCalibration[sector_csi1][ring].first << std::endl;
        // std::cout << 16*sector_csi1+ring << '\t' << Csi1RingSectorCalibration[16*sector_csi1+ring] << std::endl; 
    }
    csi1Calibration.close();

    // for(int i = 0; i<16; i++){
    //     for(int j = 0; j<16; j++){
    //         std::cout << Csi1RingSectorCalibration[i][j] << std::endl;
    //     }
    // }

    std::ifstream sd2Calibration("Sd2rInBeamCalibration.txt");
    assert(sd2Calibration.is_open());
    while(sd2Calibration >> ring >> gain >> offset) {
        Sd2rRingCalibration[ring] = std::make_pair(gain, offset);
    }
    sd2Calibration.close();

    // YY1 Sector and Ring Map
	
	for(int i = 0; i < 128; i++) {
		int sector = i/16;
		int ring = i%16;
		YY1SectorRingMap[i] = std::make_pair(sector, ring);
	}

    std::ifstream YdCalibration("Yd_4815Ped.txt");
    assert(YdCalibration.is_open());
    while(YdCalibration >> ring >> gain >> offset) {
        YdChannelCalibration[ring] = std::make_pair(gain, offset);
    }
    YdCalibration.close();

    std::ifstream YuCalibration("Yu5225_pedestal_loss_AllPar.txt");
    assert(YuCalibration.is_open());
    while(YuCalibration >> ring >> gain >> offset) {
        YuChannelCalibration[ring] = std::make_pair(gain, offset);
        // std::cout << ring << std::endl;
    }
    YuCalibration.close();
}

Bool_t Analysis::Process(Long64_t entry) {
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // When processing keyed objects with PROOF, the object is already loaded
    // and is available via the fObject pointer.
    //
    // This function should contain the \"body\" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.

    fReader.SetLocalEntry(entry);
    if(ev_num%50000 == 0) std::cout << "Current event = " << ev_num << "\r"<< std::flush;
    ev_num++;

    //Defining a structure YuDetector
    std::vector<YY1Det> YuDetector;
    for(size_t i = 0; i < *YuMul; i++) {
        if(YuChannel[i] == 82 || YuChannel[i] == 96 || YuChannel[i] == 106 || YuChannel[i] == 111) continue;
        double energy = YuEnergyRaw[i]*YuChannelCalibration[YuChannel[i]].first + YuChannelCalibration[YuChannel[i]].second;
        YY1Det hit = {*YuMul, YuChannel[i], YuRing[i], YuSector[i], YuEnergyRaw[i], energy/1000.};
        // std::cout << energy << std::endl;
        YuDetector.push_back(hit);
    }

    //Defining a structure YdDetector
    std::vector<YY1Det> YdDetector;
    for(size_t i = 0; i < *YdMul; i++) {
        if(YdChannel[i] == 0 || YdChannel[i] == 20 || YdChannel[i] == 55 || YdChannel[i] == 127) continue;
        double energy = YdEnergyRaw[i]*YdChannelCalibration[YdChannel[i]].first + YdChannelCalibration[YdChannel[i]].second;
        YY1Det hit = {*YdMul, YdChannel[i], YdRing[i], YdSector[i], YdEnergyRaw[i], energy/1000.};
        // std::cout << energy << std::endl;
        YdDetector.push_back(hit);
    }

    // std::cout << "Error before defining the structure" << std::endl;

    //Defining a structure CsI1 Detector
    std::vector<CsIDet> CsI1Detector;
    for(size_t i = 0; i < *CsI1Mul; i++){
        if(*CsI1Mul==*YdMul){
            // int csi1_chanNo = (16*CsI1Channel[i])+YdChannel[i];
            // double energy = CsI1EnergyRaw[i]*Csi1RingSectorCalibration[csi1_chanNo].first + Csi1RingSectorCalibration[csi1_chanNo].second;
            double energy = CsI1EnergyRaw[i]*Csi1RingSectorCalibration[CsI1Channel[i]][YdRing[i]].first + Csi1RingSectorCalibration[CsI1Channel[i]][YdRing[i]].second;
            // std::cout << CsI1Channel[i] << "\t" << YdRing[i] << "\t" << CsI1EnergyRaw[i] << "\t" << Csi1RingSectorCalibration[CsI1Channel[i]][YdRing[i]].first << "\t" << Csi1RingSectorCalibration[CsI1Channel[i]][YdRing[i]].second << "\t" << energy << std::endl;
            CsIDet hit = {*CsI1Mul, CsI1Channel[i], CsI1EnergyRaw[i], energy};
            CsI1Detector.push_back(hit);
        }
    }

    //Defining a structure CsI2 Detector
    std::vector<CsIDet> CsI2Detector;
    for(size_t i = 0; i < *CsI2Mul; i++){
        CsIDet hit = {*CsI2Mul, CsI2Channel[i], CsI2EnergyRaw[i]};
        CsI2Detector.push_back(hit);
    }


    //Defining a structure Sd1rDetector
    std::vector<S3Det> Sd1rDetector;
    for(size_t i = 0; i < *Sd1rMul; i++) {
        double energy = Sd1rEnergyRaw[i]*Sd1rRingCalibration[Sd1rChannel[i]].first + Sd1rRingCalibration[Sd1rChannel[i]].second;
        S3Det hit = {*Sd1rMul, Sd1rChannel[i], Sd1rEnergyRaw[i], energy/1000.};
        Sd1rDetector.push_back(hit);
    }

    //Defining a structure Sd2rDetector
    std::vector<S3Det> Sd2rDetector;
    for(size_t i = 0; i < *Sd2rMul; i++) {
        double energy = Sd2rEnergyRaw[i]*Sd2rRingCalibration[Sd2rChannel[i]].first + Sd2rRingCalibration[Sd2rChannel[i]].second;
        // std::cout << energy << '\t' << Sd2rRingCalibration[Sd2rChannel[i]].first << '\t' << Sd2rRingCalibration[Sd2rChannel[i]].second << std::endl;
        S3Det hit = {*Sd2rMul, Sd2rChannel[i], Sd2rEnergyRaw[i], energy/1000.};
        Sd2rDetector.push_back(hit);
    }



    //Sorting Yu
    if(!YuDetector.empty()) std::sort(YuDetector.begin(), YuDetector.end(), sortByEnergyYY1());

    //Sorting CsI1
    if(!CsI1Detector.empty() && !YdDetector.empty()) std::sort(CsI1Detector.begin(), CsI1Detector.end(), sortByEnergyCsI());

    //Sorting Yd
    if(!YdDetector.empty()) std::sort(YdDetector.begin(), YdDetector.end(), sortByEnergyYY1());

    //Sorting Sd1r
    if(!Sd1rDetector.empty()) std::sort(Sd1rDetector.begin(), Sd1rDetector.end(), sortByEnergyS3());

    //Sorting Sd2r
    if(!Sd2rDetector.empty()) std::sort(Sd2rDetector.begin(), Sd2rDetector.end(), sortByEnergyS3());

    // if(!pidcutDeuteron->IsInside(CsI1Detector[0].energyRaw,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))) kContinue;

    if(YuDetector.size()>1) {
        hYuEnM->Fill(YuDetector[0].energy,YuDetector[1].energy);
        hYuEnMICPID->Fill(YuDetector[0].energy,YuDetector[1].energy);
    }

    //Sd1r vs Sd2r calibrated
    if(Sd1rDetector.size() > 0 && Sd2rDetector.size() > 0 && *Sd1rMul==1 && *Sd2rMul==1) {

        double Sd1Angle = Sd1Angles[Sd1rDetector[0].channel];
        double Sd2Angle = Sd2Angles[Sd2rDetector[0].channel];

        double deutEnergyElastic = kBeam-(Sd1rDetector[0].energy+Sd2rDetector[0].energy);

        double YdAngle = Sd1Angle*sqrt((mbeam*(Sd1rDetector[0].energy+Sd2rDetector[0].energy))/(2*amu*deutEnergyElastic));
        hYdAngleElastic->Fill (YdAngle);


        hCSd1rEnRing[Sd1rDetector[0].channel]->Fill(Sd1rDetector[0].energy);
		hCSd2rEnRing[Sd2rDetector[0].channel]->Fill(Sd2rDetector[0].energy);
        hSd1rSd2r->Fill(Sd2rDetector[0].energyRaw, Sd1rDetector[0].energyRaw);
        hCSd1rSd2r->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy); //Calibrated Sd1r vs Sd2r energy with no gates

        if(pidcutElastic->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)){
            hCSd1rAngleElastic->Fill(Sd1Angle);
            hCSd2rAngleElastic->Fill(Sd2Angle);
            hCSd1rSd2rElastic->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);

            double deutEnergyElastic = kBeam-(Sd1rDetector[0].energy+Sd2rDetector[0].energy);

            double YdAngle = Sd1Angle*sqrt((mbeam*(Sd1rDetector[0].energy+Sd2rDetector[0].energy))/(2*amu*deutEnergyElastic));
            hYdAngleElastic->Fill (YdAngle);
        }

        if(*ICEnergyRaw > ICLowerBound && *ICEnergyRaw < ICUpperBound) {
            hCSd1rSd2rIC->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);  //Calibrated Sd1r vs Sd2r energy with IC gate
            if(TMath::Abs(Sd2rDetector[0].channel-Sd1rDetector[0].channel) < 2) {
                hCSd1rSd2rICCut->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);
            } //Calibrated Sd1r vs Sd2r energy with IC gate and ring vs ring correlation
            if(*YdMul==1 && YdEnergy[0] > 0) {
                hCSd1rSd2rYdIC->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);//Calibrated Sd1r vs Sd2r energy with IC gate and Yd gate
                if(TMath::Abs(Sd2rDetector[0].channel - Sd1rDetector[0].channel) < 2) {
                    hCSd1rSd2rYdICCut->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);
                } //Calibrated Sd1r vs Sd2r energy with IC and Yd gate and ring vs ring correlation
            } //end of Yd gates
        } //end of IC gates
    } //end of Sd1r vs Sd2r calibrated

    if(Sd1rDetector.size() > 0 && Sd2rDetector.size() > 0) {

        double Sd1Angle = Sd1Angles[Sd1rDetector[0].channel];
        double Sd2Angle = Sd2Angles[Sd2rDetector[0].channel];

        double AngleFactorSd1 = fabs(cos(Sd1Angle*M_PI/180.));
        double AngleFactorSd2 = fabs(cos(Sd2Angle*M_PI/180.));

        double AngleSd1_Sin = fabs(sin(Sd1Angle*M_PI/180.));

        //////////////////////////////////Reconstruction of Sd energy////////////////////////////////////////

        if(BeData){

            
            // double E_Be_Sd2_P_Target=Beam_P->AddBack(Sd2rDetector[0].energy,dSdP/AngleFactorSd2);

            // // double E_Be_Sd2_Al3_Target=Be_Al->AddBack(E_Be_Sd2_P_Target,dSdAl3);
            // // double E_Be_Sd1_Al3_Target=Be_Al->AddBack(E_Be_Sd2_Al3_Target,dSdAl3);

            // double E_Be_Sd_Al3_Target=Beam_Al->AddBack(E_Be_Sd2_P_Target, 2.*dSdAl3/AngleFactorSd2);

            // // double E_Be_Sd1_P_Target=Be_P->AddBack(E_Be_Sd1_Al3_Target,dSdP);
            // double E_Be_Sd1_P_Target=Beam_P->AddBack(E_Be_Sd_Al3_Target,dSdP/AngleFactorSd1);

            // double Sd1rE_Target=Sd1rDetector[0].energy+E_Be_Sd1_P_Target;
            // double E_Be_B_Target = Beam_B->AddBack(Sd1rE_Target,dSdB/AngleFactorSd1);
            // double E_Be_Al2_Target = Beam_Al->AddBack(E_Be_B_Target,dSdAl2/AngleFactorSd1);
            // double E_Be_SiO2_Target = Beam_SiO2->AddBack(E_Be_Al2_Target,dSdSiO2/AngleFactorSd1);
            // double E_Be_Al1_Target = Beam_Al->AddBack(E_Be_SiO2_Target,dSdAl1/AngleFactorSd1);
            // double E_Be_AgHalf = Beam_Ag->AddBack(E_Be_Al1_Target,(dAg/2)/AngleFactorSd1);//Energy loss through half of the silver

            // double E_Be_Ag_Center = E_Be_AgHalf/(AScat*pow(AngleFactorSd1+sqrt(BScat-pow(AngleSd1_Sin,2)),2));//Scattering

            // double E_Be_Ag = Beam_Ag->AddBack(E_Be_Ag_Center,(dAg/2)/AngleFactorSd1);

            // E_Be_D2_Center = Beam_D2->AddBack(E_Be_Ag,((TThickness/2)/1000.)/AngleFactorSd1);
          
            // //Energy after target and before silver

            // hSdReconAn->Fill(Sd1Angle,E_Be_Ag);
            // hSdReconT->Fill(E_Be_Ag);

            // //Energy after silver and before Sd1

            // hSdReconAfterSilverAn->Fill(Sd1Angle,E_Be_Al1_Target);
            // hSdReconAfterSilver->Fill(E_Be_Al1_Target);

            // //Energy recorded in Sd1r if there was no deadlayer 1
            // double Sd1r_EnDetect_NoDL1 = E_Be_Al1_Target-Beam_Si->CalcRemainder(E_Be_Al1_Target,dSd1Si/AngleFactorSd1);
            // hSd1NoDL1An->Fill(Sd1Angle,Sd1r_EnDetect_NoDL1);
            // hSd1NoDL1->Fill(Sd1r_EnDetect_NoDL1);


            // //Energy after Sd1 active layer and before Sd1 Dead layer 2
            
            // hSd1AfterActiveLayerAn->Fill(Sd1Angle,E_Be_Sd1_P_Target);
            // hSd1AfterActiveLayer->Fill(E_Be_Sd1_P_Target);
            

            /////////////////////Reconstruction from Sd1r Energy//////////////////////

            /*

            double Sd1rE_Be=SD1InitialEnergy(Sd1rDetector[0].energy);
            double Sd1rE_Be_B = Be_B->AddBack(Sd1rE_Be,dSdB/AngleFactorSd1);
            double Sd1rE_Be_Al2 = Be_Al->AddBack(Sd1rE_Be_B,dSdAl2/AngleFactorSd1);
            double Sd1rE_Be_SiO2 = Be_SiO2->AddBack(Sd1rE_Be_Al2,dSdSiO2/AngleFactorSd1);
            double Sd1rE_Be_Al1 = Be_Al->AddBack(Sd1rE_Be_SiO2,dSdAl1/AngleFactorSd1);
            double Sd1rE_Be_Ag = Be_Ag->AddBack(Sd1rE_Be_Al1,dAg/AngleFactorSd1);

            hSdReconSd1An->Fill(Sd1Angle,Sd1rE_Be_Ag);
            hSdReconSd1->Fill(Sd1rE_Be_Ag);

            */

            /////////////////////////////////////////////////////////////////////////

        }else{
            ///////////////////////////////12C data///////////////////////////////////
            
            // double E_C_Sd2_P_Target=Beam_P->AddBack(Sd2rDetector[0].energy,dSdP/AngleFactorSd2);

            // // double E_C_Sd2_Al3_Target=Beam_Al->AddBack(E_C_Sd2_P_Target,dSdAl3);
            // // double E_C_Sd1_Al3_Target=Beam_Al->AddBack(E_C_Sd2_Al3_Target,dSdAl3);

            // double E_C_Sd_Al3_Target=Beam_Al->AddBack(E_C_Sd2_P_Target, 2.*dSdAl3/AngleFactorSd2);

            // // double E_C_Sd1_P_Target=Beam_P->AddBack(E_C_Sd1_Al3_Target,dSdP);
            // double E_C_Sd1_P_Target=Beam_P->AddBack(E_C_Sd_Al3_Target,dSdP/AngleFactorSd1);

            // double Sd1rE_Target=Sd1rDetector[0].energy+E_C_Sd1_P_Target;
            // double E_C_B_Target = Beam_B->AddBack(Sd1rE_Target,dSdB/AngleFactorSd1);
            // double E_C_Al2_Target = Beam_Al->AddBack(E_C_B_Target,dSdAl2/AngleFactorSd1);
            // double E_C_SiO2_Target = Beam_SiO2->AddBack(E_C_Al2_Target,dSdSiO2/AngleFactorSd1);
            // double E_C_Al1_Target = Beam_Al->AddBack(E_C_SiO2_Target,dSdAl1/AngleFactorSd1);
            // double E_C_AgHalf = C_Ag->AddBack(E_C_Al1_Target,(dAg/2)/AngleFactorSd1);

            // double E_C_Ag_Center = E_C_AgHalf/(AScat*pow(AngleFactorSd1+sqrt(BScat-pow(AngleSd1_Sin,2)),2));//Scattering

            // double E_C_Ag = C_Ag->AddBack(E_C_Ag_Center,dAg/2);

            // //Energy after target and before silver

            // hSdReconAn->Fill(Sd1Angle,E_C_Ag);
            // hSdReconT->Fill(E_C_Ag);

            // //Energy after silver and before Sd1

            // hSdReconAfterSilverAn->Fill(Sd1Angle,E_C_Al1_Target);
            // hSdReconAfterSilver->Fill(E_C_Al1_Target);


            // //Energy recorded in Sd1r if there was no deadlayer 1
            // double Sd1r_EnDetect_NoDL1 = E_C_Al1_Target-Beam_Si->CalcRemainder(E_C_Al1_Target,dSd1Si/AngleFactorSd1);
            // hSd1NoDL1An->Fill(Sd1Angle,Sd1r_EnDetect_NoDL1);
            // hSd1NoDL1->Fill(Sd1r_EnDetect_NoDL1);


            // //Energy after Sd1 active layer and before Sd1 Dead layer 2
            
            // hSd1AfterActiveLayerAn->Fill(Sd1Angle,E_C_Sd1_P_Target);
            // hSd1AfterActiveLayer->Fill(E_C_Sd1_P_Target);
            

            /////////////////////Reconstruction from Sd1r Energy//////////////////////

            /*

            double Sd1rE_C=SD1InitialEnergy(Sd1rDetector[0].energy);
            double Sd1rE_C_B = Beam_B->AddBack(Sd1rE_C,dSdB/AngleFactorSd1);
            double Sd1rE_C_Al2 = Beam_Al->AddBack(Sd1rE_C_B,dSdAl2/AngleFactorSd1);
            double Sd1rE_C_SiO2 = Beam_SiO2->AddBack(Sd1rE_C_Al2,dSdSiO2/AngleFactorSd1);
            double Sd1rE_C_Al1 = Beam_Al->AddBack(Sd1rE_C_SiO2,dSdAl1/AngleFactorSd1);
            double Sd1rE_C_Ag = C_Ag->AddBack(Sd1rE_C_Al1,dAg/AngleFactorSd1);

            hSdReconSd1An->Fill(Sd1Angle,Sd1rE_C_Ag);
            hSdReconSd1->Fill(Sd1rE_C_Ag);

            */

            /////////////////////////////////////////////////////////////////////////
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////

        // std::cout << Sd2rDetector[0].energy << '\t' << E_Be_Sd2_P_Target << '\t' << E_Be_Sd_Al3_Target << '\t' << E_Be_Sd1_P_Target << std::endl;
        // std::cout << Sd1rE_Target << '\t' << E_Be_B_Target << '\t' << E_Be_Al2_Target << '\t' << E_Be_SiO2_Target << '\t' << E_Be_Al1_Target << std::endl;

        
        hSd1An->Fill(Sd1Angle,Sd1rDetector[0].energy);
        hSd2An->Fill(Sd2Angle,Sd2rDetector[0].energy);
        

        hCSd1rEn->Fill(Sd1rDetector[0].energy);
        hCSd2rEn->Fill(Sd2rDetector[0].energy);
        //std::cout << Sd1Angle << "    " << Sd2Angle << std::endl;

        if(Sd1rDetector[0].channel < 3) {
            hCSd1r_0_1_2_En->Fill(Sd1rDetector[0].energy);
        }

        if(Sd2rDetector[0].channel < 3) {
            hCSd2r_0_1_2_En->Fill(Sd2rDetector[0].energy);
        }
        if(YuDetector.empty() && YdEnergyRaw.GetSize() == 0 && CsI2EnergyRaw.GetSize() == 0) {
            //if(Sd1rDetector[0].channel==0 && (Sd2rDetector[0].channel==0 || Sd2rDetector[0].channel==1)) {
                hCSdr->Fill(Sd1rDetector[0].energy+Sd2rDetector[0].energy);
            //}
        }
    }



    ////////////////////////////////////  YU DETECTOR   ///////////////////////////////////////
    if (!YuDetector.empty() && YuDetector[0].energy > SiThresholdYu){
        //calculate the angles for the Yu detector
        yuM = Ydr0+ ( YuDetector[0].ring*YChWidth )+YChMiddle;
        TVector3 YuposM ( 0,yuM,Yuz);

        if (!random){
            YuthetaM=DetAngle[YuDetector[0].sector][YuDetector[0].ring];
        }else if (random && !BeamOffset){
            YuthetaM=ran->Uniform(Yubins[15-YuDetector[0].ring],Yubins[15-YuDetector[0].ring+1]);
        }else if (random && BeamOffset){            
            //Randomization with beam offset
            double phiRand = ran->Uniform(phi[YuDetector[0].sector]-21,phi[YuDetector[0].sector]+21);
            double rRand = ran->Uniform(r[YuDetector[0].ring]-4.9375/2, r[YuDetector[0].ring]+4.9375/2);
            double x_yuRand = rRand*sin(phiRand);
            double y_yuRand = rRand*cos(phiRand);
            TVector3 vec1Rand(0, 0, 0-Yuz);
            TVector3 vec2Rand(x_yuRand - x_target, y_yuRand - y_target, Yuz);
            double angleRand = vec1Rand.Angle(vec2Rand);
            YuthetaM=angleRand*180./M_PI;
        }

        
        YuEnergyShift=YuDetector[0].energy-shift;
            
        double YuEnergyLoss = protonELB->AddBack(YuEnergyShift, 5e-5/fabs(cos(YuthetaM*M_PI/180.)));
        YuEnergyLoss = protonELA->AddBack(YuEnergyLoss, 1e-4/fabs(cos(YuthetaM*M_PI/180.)));
        YuEnergyLoss = protonELD->AddBack(YuEnergyLoss, (TThickness/2)/1000./fabs(cos(YuthetaM*M_PI/180.)));//


        hYuEn->Fill(YuEnergyLoss);
        hYuAn->Fill(YuthetaM,YuEnergyLoss);
        Qval = ( 1+mejec/mrecoil ) * ( YuEnergyLoss) - ( 1 - mbeam/mrecoil ) * ( kBeam ) - 2 * TMath::Sqrt ( mbeam*mejec* ( YuEnergyLoss) * ( kBeam ) ) / mrecoil * TMath::Cos ( YuthetaM * TMath::Pi() / 180. );
        hQval->Fill ( Qval );
        if((Sd1rDetector.size()>0 && Sd2rDetector.size()>0) && (*ICEnergyRaw > ICLowerBound && *ICEnergyRaw < ICUpperBound)){
            if(pidcut->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)){
                hQval_cut->Fill(Qval);
                hYuEnPID->Fill(YuEnergyLoss);
                hYuAnPID_SiEn->Fill(YuthetaM,YuDetector[0].energy);
                hYuAnPID->Fill(YuthetaM,YuEnergyLoss);
            }
        }

        hQvalAn->Fill(YuthetaM,Qval);
    }

    ////////////////////////////////////  YD DETECTOR   ///////////////////////////////////////


    if (!YdDetector.empty() && YdDetector[0].energy>0.02){
        ydM = Ydr0+ ( YdDetector[0].ring *YChWidth )+YChMiddle;
		TVector3 YdposM ( 0,ydM,Ydz );
		YdposM.SetPhi ( TMath::Pi() /2-YdDetector[0].sector *TMath::Pi() /4 ); //Pi/2 because the center of the sector is at 90 degrees, Pi/4 is because there is 8 sectors so it is 2Pi/8
		YdthetaM = beam.Angle ( YdposM ) *TMath::RadToDeg();

        hYdAn->Fill ( YdthetaM,YdDetector[0].energy);
        

        if((Sd1rDetector.size()>0 && Sd2rDetector.size()>0)){
            // double Sd1Angle = Sd1Angles[Sd1rDetector[0].channel];

            if(pidcutElastic->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)){

                

                hYdAnPIDElastic->Fill ( YdthetaM,YdDetector[0].energy);
                hYdElastic[YdDetector[0].ring]->Fill (YdDetector[0].energy);
            }

            if(pidcutHorizontal1->IsInside(YdthetaM,YdDetector[0].energy)){
                hCSd1rSd2rHorizontal1->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);
            }
            if(pidcutHorizontal2->IsInside(YdthetaM,YdDetector[0].energy)){
                hCSd1rSd2rHorizontal2->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);
            }
            if(pidcutHorizontal3->IsInside(YdthetaM,YdDetector[0].energy)){
                hCSd1rSd2rHorizontal3->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy);
            }
        }


        

        if(!CsI1Detector.empty() && *CsI1Mul==*YdMul && *CsI1Mul==1){
            hYdCsI1->Fill(CsI1Detector[0].energy, (YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ));
            hYdCsI1_Tot->Fill(CsI1Detector[0].energy+YdDetector[0].energy,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ));
            hYdCsI2->Fill(CsI2Detector[0].energy, YdEnergy[0]);

            

            if(pidcutProton->IsInside(CsI1Detector[0].energyRaw,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))){
                hYdAnPIDProton->Fill ( YdthetaM,YdDetector[0].energy);
            }

            if(pidcutDeuteron->IsInside(CsI1Detector[0].energy,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))){
                hYdAnPIDDeuteron->Fill ( YdthetaM,YdDetector[0].energy);
                hYdCsI1Tot_AnPIDDeuteron->Fill ( YdthetaM,YdDetector[0].energy+CsI1Detector[0].energy);
                hCsIAnPIDDeuteron->Fill ( YdthetaM,CsI1Detector[0].energy);
            }

            if(pidcutTriton->IsInside(CsI1Detector[0].energyRaw,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))){
                hYdAnPIDTriton->Fill ( YdthetaM,YdDetector[0].energy);
            }

            if((Sd1rDetector.size()>0 && Sd2rDetector.size()>0)){
                if(pidcutProton->IsInside(CsI1Detector[0].energyRaw,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))){
                    hCSd1rSd2rProton_Yd->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy); //Calibrated Sd1r vs Sd2r energy with deuteron gate
                }

                if(pidcutDeuteron->IsInside(CsI1Detector[0].energy,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))){
                    hCSd1rSd2rDeuteron_Yd->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy); //Calibrated Sd1r vs Sd2r energy with deuteron gate
                }

                if(pidcutTriton->IsInside(CsI1Detector[0].energyRaw,(YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ))){
                    hCSd1rSd2rTriton_Yd->Fill(Sd2rDetector[0].energy, Sd1rDetector[0].energy); //Calibrated Sd1r vs Sd2r energy with deuteron gate
                }


                if(pidcutElastic->IsInside(Sd2rDetector[0].energy,Sd1rDetector[0].energy)){
                    hYdCsI1Elastic->Fill(CsI1Detector[0].energy, (YdDetector[0].energy)*TMath::Cos ( YdthetaM * TMath::Pi() / 180. ));
                    
                }
            }

            //Filling the 1D histogram for the elastic scattering of 12Be with the deuterons

            
            
        }
    } 
    

    return kTRUE;
}

void Analysis::Terminate() {
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    f_out->cd();

    //write the histograms in the output root file
    //Yu
    hSd1rSd2r->Write();
    hCSd1rSd2r->Write();
    hCSd1rAngleElastic->Write();
    hCSd2rAngleElastic->Write();
    hYdAngleElastic->Write();
    hCSd1rSd2rElastic->Write();
    hCSd1rSd2rProton_Yd->Write(); //with Yd proton cut
    hCSd1rSd2rDeuteron_Yd->Write(); //with Yd deuteron cut
    hCSd1rSd2rTriton_Yd->Write(); //with Yd triton cut
    hCSd1rSd2rHorizontal1->Write();
    hCSd1rSd2rHorizontal2->Write();
    hCSd1rSd2rHorizontal3->Write();
    hCSd1rSd2rIC->Write(); //calibrated ring vs ring with IC cut
    hCSd1rSd2rICCut->Write(); //calibrated ring vs ring with IC cut and diagonal taken into account
    hCSd1rSd2rYdIC->Write();
    hCSd1rSd2rYdICCut->Write();
    hYuAn->Write();
    hYuAnT->Write();
    hYuAn1->Write();
    hYuAnIC->Write();
    hYuAnICT->Write();
    hYuAnIC1->Write();
    hYuAnPID_SiEn->Write();
    hYuAnPID->Write();
    hYuAnPIDRad->Write();
    hYuAnPIDT->Write();
    hYuAnPID1->Write();
    hYuAnPIDYd->Write();
    hYuAnPID0_90->Write();
    hYuAnPID90_180->Write();
    hYuAnPID180_270->Write();
    hYuAnPID270_360->Write();
    hYuAnPIDRs0_7->Write();
    hYuAnPIDRs8_15->Write();


    hYdCsI1->Write();
    hYdCsI2->Write();
    hYdAn->Write();
    hYdAnPIDProton->Write();
    hYdAnPIDDeuteron->Write();
    hYdAnPIDTriton->Write();
    hYdCsI1Elastic->Write();
    hYdAnPIDElastic->Write();
    for(int i = 0; i < 16; i++) {hYdElastic[i]->Write();}
    hYdCsI1_Tot->Write();
    hCsIAnPIDDeuteron->Write();
    hYdCsI1Tot_AnPIDDeuteron->Write();

    hYuEn->Write();
    hYuEnT->Write();
    hYuEnIC->Write();
    hYuEnIC1->Write();
    hYuMul->Write();
    hYuEnPID->Write();
    hYuEn1->Write();
    hEx->Write();
    hQval->Write();
    hQval_cut->Write();
    hqvalvsEx->Write();
    hQval0_7->Write();
    hQval8_15->Write();
    hQvalT->Write();
    hQval1->Write();
    hQvalYd->Write();
    hQvalAn->Write();
    hQvalAnT->Write();
    hQval0_90->Write();
    hQval90_180->Write();
    hQval180_270->Write();
    hQval270_360->Write();
    hYuEnM->Write();
    hYuEnMICPID->Write();
    for(int i = 0; i < 8;  i++) {hYuAnPIDSec[i]->Write();}
    for(int i = 0; i < 16; i++) {hQvalR[i]->Write();}
    for(int i = 0; i < 8;  i++) {hQval2R[i]->Write();}
    for(int i = 0; i < 5;  i++) {hQval3R[i]->Write();}
    for(int i = 0; i < 4;  i++) {hQval4R[i]->Write();}
    for(int i = 0; i < 8;  i++) {hQvalS[i]->Write();}

    hAngle_1->Write();
    hAngle_2->Write();
    hAngle_3->Write();
    hAngle_4->Write();

    hSd1An->Write();
    hSd2An->Write();
    hSdReconAn->Write();
    hSdReconT->Write();

    hSdReconAfterSilverAn->Write();
    hSdReconAfterSilver->Write();
    hSd1NoDL1An->Write();
    hSd1NoDL1->Write();
    hSd1AfterActiveLayerAn->Write();
    hSd1AfterActiveLayer->Write();

    hSdReconSd1An->Write();
    hSdReconSd1->Write();
    hCSd1rEn->Write();
    hCSd2rEn->Write();
    for(int i=0; i<24; i++){hCSd1rEnRing[i]->Write();}
	for(int i=0; i<24; i++){hCSd2rEnRing[i]->Write();}

    hCSd1r_0_1_2_En->Write();
    hCSd2r_0_1_2_En->Write();

    hCSdr->Write();

    f_out->Close();

    

    // g_out->cd();

    // hYuAnPID_SiEn->Write();

    // g_out->Close();

    // delete hCSd1rSd2r;
    // delete hCSd1rSd2rIC;
    // delete hCSd1rSd2rICCut;
    // delete hCSd1rSd2rYdIC;
    // delete hCSd1rSd2rYdICCut;

    // delete hYuMul;

    // delete hYuEn;
    // delete hYuEnT;
    // delete hYuEn1;

    // delete hYuEnPID;

    // delete hYuEnIC;
    // delete hYuEnIC1;

    // delete hYuAn;
    // delete hYuAnT;
    // delete hYuAn1;

    // delete hYuAnIC;
    // delete hYuAnICT;
    // delete hYuAnIC1;

    // delete hYuAnPID;
    // delete hYuAnPIDRad;
    // delete hYuAnPIDT;
    // delete hYuAnPID1;
    // delete hYuAnPIDYd;
    // delete hYuAnPID0_90;
    // delete hYuAnPID90_180;
    // delete hYuAnPID180_270;
    // delete hYuAnPID270_360;
    // delete hYuAnPIDRs0_7;
    // delete hYuAnPIDRs8_15;
    // delete hSd1An;
    // delete hSd2An;
    // delete hSdReconAn;

    // delete hYdCsI;

    // delete hEx;
    // delete hqvalvsEx;
    // delete hQval;
    // delete hQval0_90;
    // delete hQval90_180;
    // delete hQval180_270;
    // delete hQval270_360;
    // delete hQval0_7;
    // delete hQval8_15;
    // delete hQvalT;
    // delete hQval1;
    // delete hQvalAn;
    // delete hQvalAnT;
    // delete hQvalYd;

    // delete hYuEnM;
    // delete hYuEnMICPID;

    // delete[] hYuAnPIDSec;

    // delete hCSd1rEn;
    // delete hCSd1r_0_1_2_En;

    // delete hCSd1sEn;

    // delete hCSd2rEn;
    // delete hCSd2r_0_1_2_En;

    // delete hCSd2sEn;

    // delete hCSdr;

}

double Analysis::SD1InitialEnergy(double initialEnergy) {
    double EL = 0.;
    double EH = 200.;

    double EA = (EL+EH)/2;

    double EDet = initialEnergy;

    double tolerance = 1.e-2;

    double error = 0.5;

    double ERem;
    double ER;

    long iter = 0;

    while(error > tolerance){
        if(BeData){
            ERem = Beam_Si->CalcRemainder(EA,dSd1Si);
        }else{
            ERem = Beam_Si->CalcRemainder(EA,dSd1Si);
        }
        ER = EA-ERem;
        error = fabs(EDet-ER);
        if(EDet>ER){
            EH = EA;
            EA = (EL+EH)/2;
        }else if(EDet<ER){
            EL = EA;
            EA = (EL+EH)/2;
        }
        
        iter++;
        if(iter > 1000) {
            EA = 0.;
            break;
        }
        //cout << ER << endl;
    }

    return EA;
}

std::vector<YY1Det> CheckChargeSharing(std::vector<YY1Det> detect) {
    std::vector<YY1Det> newDetector;

    int sectorHits[8];
    for(int i = 0; i < 8; i++) {
        sectorHits[i] = 0;
    }
    for(auto det : detect) {
        sectorHits[det.sector]++;
    }

    for(int i = 0; i < 8; i++) {
        if(sectorHits[i] > 1) printf("Sector %d has %d hits\n", i, sectorHits[i]);
    }

    return newDetector;
}
