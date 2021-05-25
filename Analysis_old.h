#ifndef Analysis_h
#define Analysis_h

#include <TChain.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

#include "EnergyLoss.h"

typedef struct YY1Det{
    int mult;
    int channel;
    int ring;
    int sector;
    int energyRaw;
    double energy;
} YY1Det;

typedef struct CsIDet{
    int mult;
    int channel;
    float energyRaw;
    // double energy;
} CsIDet;

typedef struct S3Det{
    int mult;
    int channel;
    int energyRaw;
    double energy;
} S3Det;

struct sortByEnergyYY1 {
    inline bool operator() (const YY1Det& EnYY1_1,
                            const YY1Det& EnYY1_2){
        return (EnYY1_1.energy>EnYY1_2.energy);
    }
};

struct sortByEnergyS3 {
    inline bool operator() (const S3Det& EnS3_1,
                            const S3Det& EnS3_2){
        return (EnS3_1.energy>EnS3_2.energy);
    }
};

class Analysis : public TSelector {
public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

    // Readers to access the data (delete the ones you do not need).
    TTreeReaderValue<Int_t> YdMul = {fReader, "YdMul"};
    TTreeReaderArray<Int_t> YdChannel = {fReader, "YdChannel"};
    TTreeReaderArray<Int_t> YdEnergyRaw = {fReader, "YdEnergyRaw"};
    TTreeReaderArray<Float_t> YdEnergy = {fReader, "YdEnergy"};
    TTreeReaderArray<Int_t> YdRing = {fReader, "YdRing"};
    TTreeReaderArray<Int_t> YdSector = {fReader, "YdSector"};
    TTreeReaderValue<Int_t> YuMul = {fReader, "YuMul"};
    TTreeReaderArray<Int_t> YuChannel = {fReader, "YuChannel"};
    TTreeReaderArray<Int_t> YuEnergyRaw = {fReader, "YuEnergyRaw"};
    TTreeReaderArray<Float_t> YuEnergy = {fReader, "YuEnergy"};
    TTreeReaderArray<Int_t> YuRing = {fReader, "YuRing"};
    TTreeReaderArray<Int_t> YuSector = {fReader, "YuSector"};
    TTreeReaderValue<Int_t> CsI1Mul = {fReader, "CsI1Mul"};
    TTreeReaderArray<Int_t> CsI1Channel = {fReader, "CsI1Channel"};
    TTreeReaderArray<Float_t> CsI1EnergyRaw = {fReader, "CsI1EnergyRaw"};
    TTreeReaderValue<Int_t> CsI2Mul = {fReader, "CsI2Mul"};
    TTreeReaderArray<Int_t> CsI2Channel = {fReader, "CsI2Channel"};
    TTreeReaderArray<Float_t> CsI2EnergyRaw = {fReader, "CsI2EnergyRaw"};
    TTreeReaderValue<Int_t> ICChannel = {fReader, "ICChannel"};
    TTreeReaderValue<Float_t> ICEnergyRaw = {fReader, "ICEnergyRaw"};
    TTreeReaderValue<Int_t> Sd1rMul = {fReader, "Sd1rMul"};
    TTreeReaderArray<Int_t> Sd1rChannel = {fReader, "Sd1rChannel"};
    TTreeReaderArray<Int_t> Sd1rEnergyRaw = {fReader, "Sd1rEnergyRaw"};
    TTreeReaderArray<Float_t> Sd1rEnergy = {fReader, "Sd1rEnergy"};
    TTreeReaderValue<Int_t> Sd1sMul = {fReader, "Sd1sMul"};
    TTreeReaderArray<Int_t> Sd1sChannel = {fReader, "Sd1sChannel"};
    TTreeReaderArray<Int_t> Sd1sEnergyRaw = {fReader, "Sd1sEnergyRaw"};
    TTreeReaderArray<Float_t> Sd1sEnergy = {fReader, "Sd1sEnergy"};
    TTreeReaderValue<Int_t> Sd2rMul = {fReader, "Sd2rMul"};
    TTreeReaderArray<Int_t> Sd2rChannel = {fReader, "Sd2rChannel"};
    TTreeReaderArray<Int_t> Sd2rEnergyRaw = {fReader, "Sd2rEnergyRaw"};
    TTreeReaderArray<Float_t> Sd2rEnergy = {fReader, "Sd2rEnergy"};
    TTreeReaderValue<Int_t> Sd2sMul = {fReader, "Sd2sMul"};
    TTreeReaderArray<Int_t> Sd2sChannel = {fReader, "Sd2sChannel"};
    TTreeReaderArray<Int_t> Sd2sEnergyRaw = {fReader, "Sd2sEnergyRaw"};
    TTreeReaderArray<Float_t> Sd2sEnergy = {fReader, "Sd2sEnergy"};
    TTreeReaderValue<Int_t> SurMul = {fReader, "SurMul"};
    TTreeReaderArray<Int_t> SurChannel = {fReader, "SurChannel"};
    TTreeReaderArray<Int_t> SurEnergyRaw = {fReader, "SurEnergyRaw"};
    TTreeReaderValue<Int_t> SusMul = {fReader, "SusMul"};
    TTreeReaderArray<Int_t> SusChannel = {fReader, "SusChannel"};
    TTreeReaderArray<Int_t> SusEnergyRaw = {fReader, "SusEnergyRaw"};
    TTreeReaderValue<Int_t> SSBADC = {fReader, "SSBADC"};
    TTreeReaderValue<Float_t> SSBEnergy = {fReader, "SSBEnergy"};
    TTreeReaderValue<Int_t> ScADC = {fReader, "ScADC"};
    TTreeReaderValue<Float_t> ScEnergy = {fReader, "ScEnergy"};
    TTreeReaderValue<Int_t> TrMul = {fReader, "TrMul"};
    TTreeReaderArray<Int_t> TrADC = {fReader, "TrADC"};
    TTreeReaderArray<Float_t> TrEnergy = {fReader, "TrEnergy"};

    Analysis(TTree * /*tree*/ =0) { }
    virtual ~Analysis() { }
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    Terminate();

    double SD1InitialEnergy(double); //Added on Feb 4th to get the initial energy from Sd1 energy

    ClassDef(Analysis, 0);

    // User Defined Variables and Functions
    bool Target;
    bool BeData;
    bool CData;
    bool TT10, TT20, TT30, TT38_6;
    bool TT43;
    bool TT45;
    bool TT47;
    bool TT49_99;
    bool TT56;
    bool TT60;
    float TThickness;

    bool SiShift;
    float shift;

    double SiThresholdYu;

    int ICLowerBound, ICUpperBound;

    bool BeamOffset;
    float x_target;
    float y_target;

    bool random;

    EnergyLoss *Beam_Al, *Beam_SiO2, *Beam_B, *Beam_Si, *Beam_P, *Beam_Ag, *Beam_D2;

    EnergyLoss* protonELB;
    EnergyLoss* protonELA;
    EnergyLoss* protonELD;

    // EnergyLoss* C_Al;
    // EnergyLoss* C_SiO2;
    // EnergyLoss* C_B;
    // EnergyLoss* C_Si;
    // EnergyLoss* C_P;
    // EnergyLoss* C_Ag;
    // EnergyLoss* C_D2;

    // EnergyLoss* Be_Al;
    // EnergyLoss* Be_SiO2;
    // EnergyLoss* Be_B;
    // EnergyLoss* Be_Si;
    // EnergyLoss* Be_P;
    // EnergyLoss* Be_Ag;
    // EnergyLoss* Be_D2;

    double dAg;
    double dSdAl1;
    double dSdSiO2;
    double dSdAl2;
    double dSdB;
    double dSd1Si;
    double dSd2Si;
    double dSdP;
    double dSdAl3;

    // Ydz distance from the target, Ydr0 inner radius, Ydz outer radius, same for the Sd
    float Ydr0;
    float Ydr1;
    float Ydz;
    float Yuz;
    float Sd1z;
    float Sd2z;
    float Sdr0;
    float Sdr1;

    double E_Be_D2_Center;

    TString path;
    TString f_out_name, f_cut_name;
    TFile *f_out;
    TFile *g_out;
    TFile *f_cut;
    TCutG *pidcut;

    TFile *f_cutTriton;
    TCutG *pidcutTriton;

    TFile *f_cutDeuteron;
    TCutG *pidcutDeuteron;

    TFile *f_cutProton;
    TCutG *pidcutProton;

    TFile *f_cutElastic;
    TCutG *pidcutElastic;

    TFile *f_cutHorizontal1;
    TCutG *pidcutHorizontal1;

    TFile *f_cutHorizontal2;
    TCutG *pidcutHorizontal2;

    TFile *f_cutHorizontal3;
    TCutG *pidcutHorizontal3;

    double Yubins[17];
    double YubinsRad[17];
    double width;

    double Ydbins[17], Ydbintemp[17];

    double Sd1Bins[25];
    double Sd1Angles[24];
    double widthSd;

    double Sd2Bins[25];
    double Sd2Angles[24];

    TH2D *hSd1rSd2r;
    TH2D *hCSd1rSd2r;
    TH1D *hCSd1rAngleElastic, *hCSd2rAngleElastic, *hYdAngleElastic;
    TH2D *hCSd1rSd2rElastic;
    TH2D *hCSd1rSd2rHorizontal1, *hCSd1rSd2rHorizontal2, *hCSd1rSd2rHorizontal3;
    TH2D *hCSd1rSd2rProton_Yd;
    TH2D *hCSd1rSd2rDeuteron_Yd;
    TH2D *hCSd1rSd2rTriton_Yd;
    TH2D *hCSd1rSd2rIC;
    TH2D *hCSd1rSd2rICCut;
    TH2D *hCSd1rSd2rYdIC;
    TH2D *hCSd1rSd2rYdICCut;

    TH1D *hYuMul;

    TH1D *hYuEn;
    TH1D *hYuEnT;
    TH1D *hYuEn1;

    TH1D *hYuEnPID;

    TH1D *hYuEnIC;
    TH1D *hYuEnIC1;

    TH2D *hYuAn;
    TH2D *hYuAnT;
    TH2D *hYuAn1;

    TH2D *hYuAnIC;
    TH2D *hYuAnICT;
    TH2D *hYuAnIC1;

    TH2D *hYuAnPID_SiEn;

    TH2D *hYuAnPID;
    TH2D *hYuAnPIDRad;
    TH2D *hYuAnPIDT;
    TH2D *hYuAnPID1;
    TH2D *hYuAnPIDYd;
    TH2D *hYuAnPID0_90;
    TH2D *hYuAnPID90_180;
    TH2D *hYuAnPID180_270;
    TH2D *hYuAnPID270_360;
    TH2D *hYuAnPIDRs0_7;
    TH2D *hYuAnPIDRs8_15;
    TH2D *hSd1An;
    TH2D *hSd2An;
    TH2D *hSdReconAn;
    TH1D *hSdReconT;
    TH2D *hSdReconSd1An;
    TH1D *hSdReconSd1;
    TH2D *hSdReconAfterSilverAn;
    TH1D *hSdReconAfterSilver;
    TH2D *hSd1NoDL1An;
    TH1D *hSd1NoDL1;
    TH2D *hSd1AfterActiveLayerAn;
    TH1D *hSd1AfterActiveLayer;

    TH2D *hYdCsI1;
    TH2D *hYdCsI1Horizontal1, *hYdCsI1Horizontal2, *hYdCsI1Horizontal3;
    TH2D *hYdCsI2;
    TH2D *hYdAn;
    TH2D *hYdAnPIDProton;
    TH2D *hYdAnPIDDeuteron;
    TH2D *hYdAnPIDTriton;
    TH2D *hYdCsI1Elastic;
    TH2D *hYdAnPIDElastic;
    TH1D *hYdElastic[16];

    TH1D *hEx;
    TH2D *hqvalvsEx;
    float qval_i, qval_f;
    int qval_bins;
    TH1D *hQval_cut;
    TH1D *hQval;
    TH1D *hQval0_90;
    TH1D *hQval90_180;
    TH1D *hQval180_270;
    TH1D *hQval270_360;
    TH1D *hQval0_7;
    TH1D *hQval8_15;
    TH1D *hQvalT;
    TH1D *hQval1;
    TH2D *hQvalAn;
    TH2D *hQvalAnT;
    TH1D *hQvalYd;

    TH2D *hYuEnM;
    TH2D *hYuEnMICPID;

    TH2D *hYuAnPIDSec[8];//YuE vs Angle with an IC and PID gate sectorwise

    TH1D *hCSd1rEnRing[24];
    TH1D *hCSd2rEnRing[24];
    TH1D *hCSd1rEn;
    TH1D *hCSd1r_0_1_2_En;

    TH1D *hCSd1sEn;

    TH1D *hCSd2rEn;
    TH1D *hCSd2r_0_1_2_En;

    TH1D *hCSd2sEn;

    TH1D *hCSdr;

    TH1D *hQvalR[16];
    TH1D *hQval2R[8];
    TH1D *hQval3R[5];
    TH1D *hQval4R[4];
    TH1D *hQvalS[8];

    TH1D *hAngle_1;
    TH1D *hAngle_2;
    TH1D *hAngle_3;
    TH1D *hAngle_4;

    long ev;
    long ev_num;

    float YChWidth;
    float YChMiddle;
    float yuM;
    float ydM;
    double YuthetaM, YuthetaMRad, YdthetaM, YuthetaSolidAngle; //angle for Yu/Yd
    TVector3 beam; // beam vector

    float SdChWidth;
    float SdChMiddle;
    float SdM;
    double SdthetaM; //angle for Sd1

    float amu;
    float massEjec;
    float kBeam; //put the correct value; beam energy 112.21 at the center of the target; 112.75 at the front of the target for 60um target
    float mAg;

    float mbeam;
    float mrecoil;
    float mejec;

    double AScat;
    double BScat;

    double Ex, Qval, QvalT, Qval1, QvalYd; //Q value variable
    std::vector<double> YuAngle;

    double phi[8];
    double r[16];
    double R[8][16];
    double DetAngle[8][16];

    int seed;
    TRandom3* ran;

    float Q13Be0;
    double YuEnergyLoss, YuEnergyShift;

    std::pair<double, double> Sd1rRingCalibration[24];
    std::pair<double, double> Sd2rRingCalibration[24];

    std::pair<double, double> YdChannelCalibration[128];
    std::pair<double, double> YuChannelCalibration[128];

    std::map<int, std::pair<int, int> > YY1SectorRingMap;

    std::vector<YY1Det> CheckChargeSharing(std::vector<YY1Det> detect);
};

#endif

#ifdef Analysis_cxx
void Analysis::Init(TTree *tree) {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    fReader.SetTree(tree);
}

Bool_t Analysis::Notify() {
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


#endif // #ifdef Analysis_cxx
