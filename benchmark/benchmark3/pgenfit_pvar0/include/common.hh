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

//! \file common.hh
//! \brief Common class containing common data structures
#ifndef common_h
#define common_h 1

#include <string>
#include "TTree.h"
#include <vector>
#include <list>

#define NCPUS_UNBINFIT 24

#define EVAL_FAST
#define PATHFLOW
//set if you wish parent neutron efficiency to be distributed from 40 to 68%
#define PARENT_NEUEFF_UNIFORM

#define ENTRYLIMIT 5000 //set negative for not limiting the entries by default
#define STARTFIT 0.08

// set if production rate of isiomeric state equal 1
#define ISOMER_SUM_UNITY

// set if flat backgrounds is used
#define FLAT_BACKGROUNDS

//#define DEBUG

#define kmaxndecay 200
#define kmaxpaths 200
#define kmaxparms 150
#define kmaxnri 100

//! make path
typedef struct {
    // idendification
    Int_t id;
    Int_t z;
    Int_t n;
    TString name;

    // decay properies
    Double_t decay_hl;
    Double_t decay_lamda;

    Double_t decay_p0n;//decay to several isomerics states or ground state
    Double_t decay_p1n;
    Double_t decay_p2n;


    Double_t decay_hlerr;
    Double_t decay_lamdaerr;
    Double_t decay_p0nerr;
    Double_t decay_p1nerr;
    Double_t decay_p2nerr;

    Double_t decay_hlup;
    Double_t decay_lamdaup;
    Double_t decay_p0nup;
    Double_t decay_p1nup;
    Double_t decay_p2nup;


    Double_t decay_hllow;
    Double_t decay_lamdalow;
    Double_t decay_p0nlow;
    Double_t decay_p1nlow;
    Double_t decay_p2nlow;

    Double_t population_ratio;//isomer population ratio
    Double_t population_ratioerr;
    Double_t population_ratioup;
    Double_t population_ratiolow;

    Double_t neueff;
    Double_t neuefferr;
    Double_t neueffup;
    Double_t neuefflow;

    // fit options 0-vary, 1-fix with error propagation, 2-fix without error propagation
    Int_t is_decay_hl_fix;
    Int_t is_decay_lamda_fix;
    Int_t is_decay_p0n_fix;
    Int_t is_decay_p1n_fix;
    Int_t is_decay_p2n_fix;
    Int_t is_population_ratio_fix;
    Int_t is_neueff_fix;

    Int_t gspatner; //ground statte patner, set negative if gs itself

    // paths to this ri
    std::vector< std::vector<Int_t> > path;
    std::vector< std::vector<Int_t> > nneupath;



    // stuffs for simulation
    Int_t sim_neumult;
    Double_t sim_T;
    Bool_t sim_ispopulated;
} MemberDef;

//! reduced paths
typedef struct path_str
{
    Int_t nri;
    Int_t npaths;
    Int_t ndecay[kmaxpaths];
    Int_t decaymap[kmaxpaths][kmaxndecay];
    Int_t nneu[kmaxpaths][kmaxndecay];
    Bool_t ispathhasflow[kmaxpaths];

    //! stuffs for faster calculation
    Int_t nri5;
    Int_t nri4;
    Int_t nri3;
    Int_t nri2;

    //! stuffs for isomer
    Int_t nisomers;
    Int_t isomer_gs_index[10];
    Int_t isomer_ex_index[10];
} path;


#endif
