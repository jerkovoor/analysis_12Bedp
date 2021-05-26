#include "Analysis.h"

TChain* MakeChain();

int main() {
    TChain *chain = MakeChain();

    // Mode 0
    auto *analysis12Be = new Analysis();
    chain->Process(analysis12Be);

    return 0;
}

TChain* MakeChain() {
    auto *chain = new TChain("AutoTree");

    bool BeData = 0;
    bool Target = 1;
    bool Server = 1;


    TString InPath;

    if(Server){
        InPath = "/home/reactions/Decoded_root";
    }else{
        InPath = "/mnt/c/Users/jerom/Desktop/12Be_exp/Analysis/NewDecodeJosh";
    }


    TString PathToFilesBeT = InPath + "/Be/NewDecode_";
    TString PathToFilesBeNT = InPath + "/Be_notarget/NewDecode_";
    TString PathToFilesCT = InPath + "/C/NewDecode_";
    TString PathToFilesCNT = InPath + "/C_notarget/NewDecode_";


    /*
    std::vector<int> run_{5021};

    for(size_t i = 0; i < run_.size(); i++) {
        chain->Add(PathToFiles + Form("%d.root", run_[i]));
    }
    */


   if(BeData){
       if(Target){
           //Be target data
            for(int run_num=5021;run_num<5096;run_num++){
                    if(run_num==5026||run_num==5040||run_num==5043||run_num==5046||run_num==5047||run_num==5059||run_num==5062||run_num==5063||run_num==5093||run_num==5099||run_num==5101){
                        continue;
                    }else{
                        chain->Add(PathToFilesBeT + Form("%d.root", run_num));
                    }
            }

            for(int run_num=5096;run_num<5113;run_num++){
                    if(run_num==5093||run_num==5099||run_num==5101){
                        continue;
                    }else{
                        chain->Add(PathToFilesBeT + Form("%d.root", run_num));
                    }
            }

            for(int run_num=5115;run_num<5133;run_num++){
                chain->Add(PathToFilesBeT + Form("%d.root", run_num));
            }

            for(int run_num=5180;run_num<5223;run_num++){
                    if(run_num==5187||run_num==5200){
                        continue;
                    }else{
                        chain->Add(PathToFilesBeT + Form("%d.root", run_num));
                    }
            }
       }else{
           //Be noTarget data
            chain->Add (PathToFilesBeNT + Form("%d.root", 5018));
            chain->Add (PathToFilesBeNT + Form("%d.root", 5135));

            for(int run_num=5140;run_num<5178;run_num++){
                    if(run_num==5141||run_num==5149||run_num==5150||run_num==5153||run_num==5160||run_num==5175){
                        continue;
                    }else{
                        chain->Add(PathToFilesBeNT + Form("%d.root", run_num));
                    }
            }
       }
   }else{
       if(Target){
           //C target data
            for(int run_num=5001;run_num<5010;run_num++){
                    if(run_num==5005||run_num==5008){
                        continue;
                    }else{
                        chain->Add(PathToFilesCT + Form("%d.root", run_num));
                    }
            }
       }else{
            //C no target data
            chain->Add(PathToFilesCNT + "4992.root");
       }
   }
     

    return chain;
}
