#include <iostream>
#include <string>
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/PointSetInterval.h"
#include "RooFitResult.h"
#include "TH1F.h"
#include "RooAbsArg.h"
#include "RooWorkspace.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
using namespace RooFit;
using namespace RooStats;


void generatedata(RooWorkspace *wks){
    // Set range of observable
   Double_t low = 80, high = 200;

    // Create a variable for the observable (invariant mass)
    RooRealVar invMass("invMass", "M_{inv}", low, high, "GeV");

    invMass.setBins(50);

    //----------------------------------- Define the PDFs for background and signal

    // Create background dataset (exponential)
    RooRealVar alpha("alpha", "#alpha", -0.01,-0.2,0.01);
    RooExponential background("background", "Background PDF", invMass, alpha);

    // Create signal dataset (gaussian)
    RooRealVar mean("mean", "Signal Mean", 125, 90, 160);
    RooRealVar sigma("sigma", "Signal Sigma", 10, 0, 20);
    RooGaussian signal("signal", "Signal PDF", invMass, mean, sigma);

    mean.setConstant();
    sigma.setConstant();

    //-----------------------------------Define signal parameters for the model

    // Introduce mu: the signal strength in units of the expectation.
    // eg. mu = 1 is the SM, mu = 0 is no signal, mu=2 is 2x the SM
    RooRealVar mu("mu", "signal strength in units of SM expectation", 1, 0., 2);

    RooRealVar fsigExpected("fsigExpected", "expected fraction of signal events", .05, 0., 1);
    fsigExpected.setConstant(); // use mu as main parameter, so fix this.

    RooRealVar ratioSigEff("ratioSigEff", "ratio of signal efficiency to nominal signal efficiency", 1.0, 0., 2);
    ratioSigEff.setConstant(kTRUE);

    RooProduct fsig("fsig", "fraction of signal events", RooArgSet(mu, ratioSigEff, fsigExpected));

    //-----------------------------------Define model from data and background PDFs
    
    RooAddPdf model("model", "Data", RooArgList(signal, background), fsig);

    //----------------------------------- Generate events based on PDFs for background, signal and model respectively

    //Toy MC generation - DataSet (unbinned)
    RooDataSet *bkgData = background.generate(RooArgSet(invMass), 10000);
    RooDataSet *sigData = signal.generate(RooArgSet(invMass), 300);
    RooDataSet *Data = model.generate(RooArgSet(invMass), 5000);

    //----------------------------------- Generate binned events and HistPDFs

    //Toy MC generation - DataHist (binned)
    RooDataHist* bkgDataHist = bkgData->binnedClone() ;
    RooDataHist* sigDataHist = sigData->binnedClone() ;
    RooDataHist* DataHist = Data->binnedClone() ;

    //Hist PDF from toy MC
    RooHistPdf bkgHistPdf("bkg","bkg",invMass,*bkgDataHist,0);
    RooHistPdf sigHistPdf("sig","sig",invMass,*sigDataHist,0);

    //Create binned_model 
    RooAddPdf binned_model("binned model","binned model",RooArgList(sigHistPdf,bkgHistPdf),fsig);

    //----------------------------------- Import model and data to workspace for Hypothesis test

    wks->import(model);
    wks->import(*Data, Rename("data")); 

    //----------------------------------- Make plots for signal and background 

    //Create a ROOT Canvas
    TCanvas *bkg_can = new TCanvas();
    //Create a plot of generated data superimposed to the function for background
    RooPlot *bkg_plot = invMass.frame();
    bkgData->plotOn(bkg_plot);
    background.plotOn(bkg_plot, RooFit::LineColor(kBlue));
    //background.paramOn(bkg_plot);

    //Draw components and save as a .png
    bkg_plot->SetTitle("Background PDF superimposed to generated data");
    bkg_plot->Draw();
    bkg_can->Draw();
    bkg_can->SaveAs("background.png");

    //Create a ROOT Canvas
    TCanvas *sig_can = new TCanvas();

    //Create a plot of generated data superimposed to the function for signal
    RooPlot *sig_plot = invMass.frame();
    sigData->plotOn(sig_plot);
    signal.plotOn(sig_plot, RooFit::LineColor(kRed));
    //signal.paramOn(sig_plot);

    //Draw components and save as a .png
    sig_plot->SetTitle("Signal PDF superimposed to generated data");
    sig_plot->Draw();
    sig_can->Draw();
    sig_can->SaveAs("signal.png");

    //Create a ROOT Canvas
    TCanvas *bkg_binned_can = new TCanvas();

    //Create a plot of the generated data superimposed to the background model
    RooPlot *bkg_binned_plot = invMass.frame();
    bkgDataHist->plotOn(bkg_binned_plot);
    bkgHistPdf.plotOn(bkg_binned_plot, RooFit::DrawOption("F"), RooFit::FillColor(kBlue), RooFit::FillStyle(1001));
    bkgHistPdf.paramOn(bkg_binned_plot);

    //Draw components and save as a .png
    bkg_binned_plot->SetTitle("Background HistPDF superimposed to generated data");
    bkg_binned_plot->Draw();
    bkg_binned_can->Draw();
    bkg_binned_can->SaveAs("background_binned.png");

    //Create a ROOT Canvas
    TCanvas *sig_binned_can = new TCanvas();

    //Create a plot of the generated data superimposed to the background model
    RooPlot *sig_binned_plot = invMass.frame();
    sigDataHist->plotOn(sig_binned_plot);
    sigHistPdf.plotOn(sig_binned_plot, RooFit::DrawOption("F"), RooFit::FillColor(kRed), RooFit::FillStyle(1001));
    //sigHistPdf.paramOn(sig_binned_plot);

    //Draw components and save as a .png
    sig_binned_plot->SetTitle("Signal HistPDF superimposed to generated data");
    sig_binned_plot->Draw();
    sig_binned_can->Draw();
    sig_binned_can->SaveAs("signal_binned.png");

    //----------------------------------- Make plots for model

    binned_model.fitTo(*DataHist);

    //Plot the binned_model (I am using the "Data" toyMC instead of the newly generated "DataHist" because it seems to be the same, not sure if it is correct)
    TCanvas *model_binned_can = new TCanvas();
    RooPlot *model_binned_plot = invMass.frame();
    DataHist->plotOn(model_binned_plot,Name("data"));
    binned_model.plotOn(model_binned_plot,Name("model"), RooFit::DrawOption("F"), RooFit::FillColor(kRed), RooFit::FillStyle(1001));
    binned_model.plotOn(model_binned_plot,Name("background"), RooFit::Components(bkgHistPdf), RooFit::DrawOption("F"), RooFit::LineStyle(2) , RooFit::FillStyle(1001), RooFit::FillColor(kBlue));
    binned_model.paramOn(model_binned_plot);

    //Create Legend
    TLegend *leg = new TLegend(0.65,0.73,0.86,0.87);
    leg->SetFillColor(kWhite);
    leg->SetLineColor(kBlack);
    leg->AddEntry(model_binned_plot->findObject("data"), "Data", "P");
    leg->AddEntry(model_binned_plot->findObject("model"),"Signal + background","F");
    leg->AddEntry(model_binned_plot->findObject("background"), "Background only", "F");
    
    //Draw components and save as .png
    model_binned_plot->SetTitle("An example fit to the signal + background model for binned dataset");
    model_binned_plot->Draw();
    model_binned_can->Draw();
    leg->Draw();
    model_binned_can->SaveAs("model_binned.png");


    mu.setConstant(kFALSE);
    
    model.fitTo(*Data, Save(kTRUE), Minos(kFALSE), Hesse(kFALSE), PrintLevel(-1));
    
    // plot sig candidates, full model, and individual components
    TCanvas *sbmodel_can =new TCanvas();
    RooPlot *sbmodel_plot = invMass.frame();
    Data->plotOn(sbmodel_plot, Name("data"));
    model.plotOn(sbmodel_plot, Name("model"), LineColor(kViolet));
    model.plotOn(sbmodel_plot, Name("signal only"), Components(signal), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(sbmodel_plot, Name("background"), Components(background), LineStyle(kDashed), LineColor(kBlue));
    model.paramOn(sbmodel_plot);

    //Create Legend
    TLegend *leg1 = new TLegend(0.65,0.73,0.86,0.87);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kBlack);
    leg1->AddEntry(sbmodel_plot->findObject("data"), "Data", "P");
    leg1->AddEntry(sbmodel_plot->findObject("model"),"Signal + background","L");
    leg1->AddEntry(sbmodel_plot->findObject("background"), "Background only", "L");
    leg1->AddEntry(sbmodel_plot->findObject("signal only"), "Signal only", "L");

    //Draw components and save as .png
    sbmodel_plot->SetTitle("An example fit to the signal + background model");
    sbmodel_plot->Draw();
    sbmodel_can->Draw();
    leg1->Draw();
    sbmodel_can->SaveAs("signal+background_model.png");

    // Set signal fraction to be a constant and speciffically 0 
    mu.setVal(0);          
    mu.setConstant(kTRUE); 
    
    model.fitTo(*Data, Save(kTRUE), Minos(kFALSE), Hesse(kFALSE), PrintLevel(-1));
    
    // Plot signal candidates with background model and components
    TCanvas *bmodel_can =new TCanvas();
    RooPlot *bmodel_plot = invMass.frame();
    Data->plotOn(bmodel_plot,Name("data"), DataError(RooAbsData::SumW2));
    model.plotOn(bmodel_plot,Name("model"), LineColor(kViolet));
    model.plotOn(bmodel_plot,Name("background"), Components(background), LineStyle(kDashed), LineColor(kBlue));
    model.paramOn(bmodel_plot);

    //Create Legend
    TLegend *leg2 = new TLegend(0.65,0.73,0.86,0.87);
    leg2->SetFillColor(kWhite);
    leg2->SetLineColor(kBlack);
    leg2->AddEntry(bmodel_plot->findObject("data"), "Data", "P");
    leg2->AddEntry(bmodel_plot->findObject("model"),"Signal + background","L");
    leg2->AddEntry(bmodel_plot->findObject("background"), "Background only", "L");

    //Draw components and save as .png
    bmodel_plot->SetTitle("An example fit to the background-only model");
    bmodel_plot->Draw();
    bmodel_can->Draw();
    leg2->Draw();
    bmodel_can->SaveAs("backgroundonly.png");
}

void DoHypothesisTest(RooWorkspace *wks){
    
    // Use a RooStats ProfileLikleihoodCalculator to do the hypothesis test.
    ModelConfig model;
    model.SetWorkspace(*wks);
    model.SetPdf("model");
    
    // Load from workspace the data and the mu parameter. Set model and set mu as poi
    ProfileLikelihoodCalculator plc;
    plc.SetData(*(wks->data("data")));
    plc.SetModel(model);

    RooRealVar *mu = wks->var("mu");
    RooArgSet poi(*mu);

    // Here we explicitly set the value of the parameters for the null.
    // We want no signal contribution, mu = 0
    RooArgSet *nullParams = (RooArgSet *)poi.snapshot();
    nullParams->setRealValue("mu", 0);
    
    // Set the other parameters as nuisance
    plc.SetNullParameters(*nullParams);
    
    // Get the result of the hypothesis test from the calculator.
    HypoTestResult *htr = plc.GetHypoTest();

    //Calculate p-value and significance
    double p_value = htr->NullPValue();
    double Significance = htr->Significance();

    cout << "-------------------------------------------------" << endl;
    cout << "The p-value for the null hypothesis is " << p_value << endl;
    cout << "Which coresponds  to a significance of " << Significance << " sigma" <<endl; 
    cout << "-------------------------------------------------\n\n" << endl;

    htr->Print();




}

int Example1() {

    // Create a workspace to manage the project.
    RooWorkspace *wspace = new RooWorkspace("myWS");

    generatedata(wspace);

    DoHypothesisTest(wspace);


    return 0;
}