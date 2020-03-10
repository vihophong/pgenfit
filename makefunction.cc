#include <iostream>
#include "unbinfit.hh"

int main(int argc, char *argv[])
{
    if (argc==2) {
        unbinfit* fit=new unbinfit;
        //keep default start time (0.08s)
        char inputRootFile[1000];
        sprintf(inputRootFile,"dummy.root");
        fit->Init(argv[1],inputRootFile);
        fit->generateRoofitEvaluate();
    }else{
        std::cout<<"check inputs!"<<std::endl;
        return 0;
    }
    return 0;
}
