#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"

using namespace RooFit; 


void HiggsModel()
{ 
   //Set the number of signal and background events 
   int nsig = 100; 
   int nbkg = 1000;   

   RooWorkspace w("w"); 
   w.factory("Exponential:bkg_pdf(x[80,200], a[-0.01,-0.2,0.01])");
   w.factory("Gaussian:sig_pdf(x, mass[125], sigma[10])");

   //Create extended model 
   w.factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[0,10000]*bkg_pdf)"); 

   //Extract from workspace the pdf and the observable variable 
   RooAbsPdf * pdf = w.pdf("model");
   RooRealVar * x = w.var("x");  

   // set the desired value of signal and background events
   w.var("nsig")->setVal(nsig);
   w.var("nbkg")->setVal(nbkg);

   // generate the data

   // use fixed random numbers for reproducibility (use 0 for changing every time)
   RooRandom::randomGenerator()->SetSeed(111);

   // fix number of bins to 50 to plot or to generate data (default is 100 bins) 
   x->setBins(50);

   RooDataSet * data = pdf->generate( *x);  // will generate according to total S+B events
   //RooDataSet * data = pdf->generate( *x, AllBinned());  // will generate accordint to total S+B events
   data->SetName("data");
   w.import(*data);

   data->Print(); 

   //Create Canvas 
   TCanvas *can = new TCanvas();

   RooPlot * plot = x->frame(Title("Gaussian Signal over Exponential Background"));
   data->plotOn(plot, Name("data"));
   plot->Draw();

   RooFitResult * r = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
   r->Print();

   pdf->plotOn(plot, Name("model"), RooFit::LineColor(kViolet));
   //draw the two separate pdf's
   pdf->plotOn(plot, Name("background"), RooFit::Components("bkg_pdf"), RooFit::LineColor(kBlue) ,RooFit::LineStyle(kDashed) );
   pdf->plotOn(plot, Name("signal only"), RooFit::Components("sig_pdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );

   pdf->paramOn(plot,Layout(0.5,0.9,0.85));

   //Create Legend
   TLegend *leg = new TLegend(0.65,0.73,0.86,0.87);
   leg->SetFillColor(kWhite);
   leg->SetLineColor(kBlack);
   leg->AddEntry(plot->findObject("data"), "Data", "P");
   leg->AddEntry(plot->findObject("model"),"Signal + background","L");
   leg->AddEntry(plot->findObject("background"), "Background only", "L");
   leg->AddEntry(plot->findObject("signal only"), "Signal only", "L");

   plot->Draw();
   can->Draw();
   leg->Draw();
   can->SaveAs("Signal+background_model.png");

   //Create the ModelConfig in order to use later 
   RooStats::ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*pdf);
   mc.SetParametersOfInterest(*w.var("nsig"));
   mc.SetObservables(*w.var("x"));
   // define set of nuisance parameters
   w.defineSet("nuisParams","a,nbkg");

   mc.SetNuisanceParameters(*w.set("nuisParams"));

   // import model in the workspace 
   w.import(mc);

   // write the workspace in the file
   TString fileName = "HiggsModel.root";
   w.writeToFile(fileName,true);
   cout << "model written to file " << fileName << endl;
}
