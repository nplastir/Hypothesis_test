#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"

using namespace RooFit; 


void HiggsHistModel()
{ 
   RooWorkspace wsim("wsim"); 
   wsim.factory("Exponential:bkg_pdf(x[40,400], a[-0.01,-10,0])");
   wsim.factory("Gaussian:sig_pdf(x, mass[125, 80 , 400], sigma[5,1,10])");

   //Create binned data
   RooDataHist* hist_sig = wsim.pdf("sig_pdf")->generateBinned(*wsim.var("x"),200) ;
   RooDataHist* hist_bkg = wsim.pdf("bkg_pdf")->generateBinned(*wsim.var("x"),10000);

   //Create mock data with mu = 1.5
   wsim.factory("expr::S('mu*Snom',mu[1.5],Snom[200])") ;
   wsim.factory("SUM::model(S*sig_pdf,Bnom[10000]*bkg_pdf)") ;
   RooDataHist* hist_data = wsim.pdf("model")->generateBinned(*wsim.var("x")) ;

   //Create the binned likelihood
   RooWorkspace w("w") ;

   w.import(*hist_sig,RooFit::Rename("template_sig")) ;
   w.import(*hist_bkg,RooFit::Rename("template_bkg")) ;
   w.import(*hist_data,RooFit::Rename("observed_data")) ;

   w.factory("HistFunc::sig(x,template_sig)") ;
   w.factory("HistFunc::bkg(x,template_bkg)") ;

   w.factory("binw[0.277]") ; // bin width == 1/(400-30)
   w.factory("L[1]") ; // L is the luminosity ratio data to simulation. E.g. if L(simul) = 3x L(data), L=0.33 (so that all prediction are scaled down to what is expected for the data)
   w.factory("expr::S('mu*L*binw',mu[1,-1,6],L,binw[0.277])") ;
   w.factory("expr::B('Bscale*L*binw',Bscale[0,6],L,binw)") ;
   w.factory("ASUM::model(S*sig,B*bkg)") ;

   w.pdf("model")->fitTo(*hist_data) ;

   TCanvas* c1 = new TCanvas() ;
   RooPlot* frame = w.var("x")->frame() ;
   hist_data->plotOn(frame,Name("data")) ;
   w.pdf("model")->plotOn(frame, Name("model"), RooFit::LineColor(kViolet)) ;
   w.pdf("model")->plotOn(frame,Name("background"), RooFit::Components("bkg"), RooFit::LineColor(kBlue) ,RooFit::LineStyle(kDashed) );
   w.pdf("model")->plotOn(frame, Name("signal only"), RooFit::Components("sig"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );

   //Create Legend
   TLegend *leg = new TLegend(0.65,0.73,0.86,0.87);
   leg->SetFillColor(kWhite);
   leg->SetLineColor(kBlack);
   leg->AddEntry(frame->findObject("data"), "Data", "P");
   leg->AddEntry(frame->findObject("model"),"Signal + background","L");
   leg->AddEntry(frame->findObject("background"), "Background only", "L");
   leg->AddEntry(frame->findObject("signal only"), "Signal only", "L");

   frame->Draw();
   leg->Draw();
   c1->Draw() ;
   c1->SaveAs("Signal+background_model.png");


   //Create the ModelConfig in order to use later 
   RooStats::ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
   mc.SetParametersOfInterest(*w.var("mu"));
   mc.SetObservables(*w.var("x"));
   mc.SetNuisanceParameters(*w.var("Bscale"));
   // // define set of nuisance parameters
   // w.defineSet("nuisParams","a,nbkg");
   // mc.SetNuisanceParameters(*w.set("nuisParams"));

   w.var("mu")->setVal(1) ;
   mc.SetSnapshot(*w.var("mu"));
   mc.Print() ;

   // import model in the workspace 
   w.import(mc);

   // write the workspace in the file
   TString fileName = "HiggsHistModel.root";
   w.writeToFile(fileName,true);
   cout << "model written to file " << fileName << endl;
}
