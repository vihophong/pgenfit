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
#include "simulation.hh"
#include "TFile.h"

int main(int argc, char *argv[])
{
    if (argc==4) {
        simulation* sim=new simulation(argv[1]);

        sim->readSimulationParameters(argv[2]);
        TFile* fout=new TFile(argv[3],"recreate");
        fout->cd();
        sim->setRandomSeed(0);
        sim->BookSimulationTree();
        sim->BookCorrelationTree();
        sim->runSimulation();
        TTree* treeion=sim->getIonSimulationTree();
        TTree* treebeta=sim->getBetaSimulationTree();
        TTree* treeneutron=sim->getNeutronSimulationTree();
        TTree* treecorr=sim->getCorrelationTree();
        TTree* treemlh=sim->getMLHTree();
        sim->fillTreeData();
        sim->correlateData();
        treeion->Write();
        treebeta->Write();
        treeneutron->Write();
        treecorr->Write();
        treemlh->Write();
        sim->writeMLHHistos();
        fout->Close();
    }else{
        std::cout<<"check inputs!"<<std::endl;
        return 0;
    }
    return 0;
}
