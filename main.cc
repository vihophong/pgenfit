//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * Copyright@2019 Vi Ho Phong, email: phong@ribf.riken.jp           *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications.                    *
// ********************************************************************
//

#include <iostream>
#include <fstream>
#include "unbinfit.hh"

int main(int argc, char *argv[])
{
    int c;
    extern char *optarg;
    std::string plotrange_str;
    while ((c = getopt(argc, argv, "r:")) != EOF)
    {

        switch (c)
        {
        case 'r':   // range
            plotrange_str = std::string(optarg);
            break;

        default:
            break;
        }
    }
    double plotrangelow=-9999, plotrangehi=-9999;
    if (!plotrange_str.empty()) {
        std::replace(plotrange_str.begin(), plotrange_str.end(), ',', ' ');
        std::istringstream iss(plotrange_str);
        if (!(iss >> plotrangelow >> plotrangehi)) {
            std::cerr << "Invalid format for -c. Use x,y with float values\n";
            return 1;
        }
    }

    /*
    unbinfit* fit=new unbinfit;
    char inputRootFile[1000];
    sprintf(inputRootFile,"testdata.root");
    fit->Init(argv[1],inputRootFile);
    fit->generateRoofitEvaluate();
    */
    std::cout << "plot range low = " << plotrangelow << ", plot range high = " << plotrangehi << "\n";
    if (argc==2) {
        unbinfit* fit=new unbinfit;
        char inpparms[1000];
        char inputRootFile[1000];
        char inputeffparms[1000];
        sprintf(inpparms,"parmsex.txt");
        sprintf(inputeffparms,"effparmsex.txt");
        sprintf(inputRootFile,"testdata.root");
        //keep default start time (0.08s)
        fit->setStartTime(0.04);
        fit->setNBinHists(200);
        fit->SetRandomSeed(4357);//must be set before Init
        fit->setEntriesLimit(5000);
        fit->Init(inpparms,inputRootFile);
        fit->setInputEffParms(inputeffparms);
        fit->setOutputFile(argv[1]);
        char outputTextFile[1000];
        sprintf(outputTextFile,"%s.txt",argv[1]);
        std::ofstream ofs(outputTextFile);
        for (Int_t i=0;i<argc;i++){
            ofs<<argv[i]<<"\t";
        }
        ofs<<std::endl;
        if (plotrangelow>-9999){
            fit->setPlotTimeRange(plotrangelow,plotrangehi);
        }
        fit->Run();
    }else if(argc==9){
        unbinfit* fit=new unbinfit;
        fit->setStartTime(atof(argv[4]));
        fit->setNBinHists(atoi(argv[5]));
        fit->SetRandomSeed(atoi(argv[8]));//must be set before Init
        fit->Init(argv[1],argv[2]);
        fit->setOutputFile(argv[3]);
        fit->setNumberOfMC(atoi(argv[6]));
        fit->setInputEffParms(argv[7]);
        if (plotrangelow>-9999){
            fit->setPlotTimeRange(plotrangelow,plotrangehi);
        }
        fit->Run();
    }else if (argc==10){
        unbinfit* fit=new unbinfit;
        fit->setStartTime(atof(argv[4]));
        fit->setTimeRange(atof(argv[9]));
        fit->setNBinHists(atoi(argv[5]));
        fit->SetRandomSeed(atoi(argv[8]));//must be set before Init
        fit->Init(argv[1],argv[2]);
        fit->setOutputFile(argv[3]);
        fit->setNumberOfMC(atoi(argv[6]));
        fit->setInputEffParms(argv[7]);
        if (plotrangelow>-9999){
            fit->setPlotTimeRange(plotrangelow,plotrangehi);
        }
        fit->Run();
    }else{
        std::cout<<"check inputs!"<<std::endl;
        std::cout<<"example: ./runsinglefit.sh parms/Cd133parms.txt lowin/Cd133.root outFit/Cd133out_unbin_withMC.root 0.08 500 525 effparmsex.txt 4357"<<std::endl;
        return 0;
    }

    std::cout<<"Fitting done, output file at "<<argv[3]<<std::endl;
    return 0;
}
