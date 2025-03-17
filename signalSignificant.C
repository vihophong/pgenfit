#include <iostream>
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/LikelihoodInterval.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

void signalSignificant() {
    // Observable: Mass variable
    RooRealVar mass("mass", "Mass", 0., 160.);

    // Signal: Gaussian (e.g., Higgs peak at 125 GeV)
    RooRealVar mean("mean", "Mean of Gaussian", 125., 100., 150.);
    RooRealVar sigma("sigma", "Width of Gaussian", 2., 0.1, 10.);
    RooGaussian signal("signal", "Signal PDF", mass, mean, sigma);

    // Background: Exponential
    RooRealVar tau("tau", "Decay parameter", -0.1, -1., 0.);
    RooExponential background("background", "Background PDF", mass, tau);

    // Yields: Number of signal and background events
    RooRealVar nsig("nsig", "Signal yield", 5, 0, 1000);
    RooRealVar nbkg("nbkg", "Background yield", 500, 0, 10000);

    // Total model: Signal + Background
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // Generate Toy Data (Replace with real data if available)
    RooDataSet *data = model.generate(RooArgSet(mass), 1000);

    // Fit model to data
    model.fitTo(*data, PrintLevel(-1));

    // Print fit results
    nsig.Print();
    nbkg.Print();

    // === Compute Significance using Profile Likelihood ===
    ProfileLikelihoodCalculator plc(*data, model, RooArgSet(nsig));
    LikelihoodInterval* interval = plc.GetInterval();
    double upperLimit = interval->UpperLimit(nsig);
    cout << "Upper Limit on nsig: " << upperLimit << endl;

    // Convert Upper Limit to significance (sigma)
    double significance = upperLimit / nsig.getError();  // Approximate significance
    cout << "Signal significance (Profile Likelihood): " << significance << " sigma" << endl;
    delete data;
}

