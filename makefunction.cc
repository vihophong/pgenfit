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
