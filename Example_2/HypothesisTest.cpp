using namespace RooStats;
using namespace RooFit;

void HypothesisTest( const char* filename =  "HiggsModel.root", 
                     const char* workspaceName = "w",
                     const char* modelConfigName = "ModelConfig",
                     const char* dataName = "data" )
{
    /////////////////////////////////////////////////////////////
    // First part is just to access the workspace file 
    ////////////////////////////////////////////////////////////

    // open input file 
    TFile *file = TFile::Open(filename);
    if (!file) return;

    // get the workspace out of the file
    RooWorkspace* w = (RooWorkspace*) file->Get(workspaceName);


    // get the data  out of the file
    RooAbsData* data = w->data(dataName);

    // Get the ModelConfig out of the file
    ModelConfig*  sbModel = (RooStats::ModelConfig*) w->obj(modelConfigName);
    sbModel->SetName("S+B Model");      
    RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
    poi->setVal(50);  // set POI snapshot in S+B model for expected significance
    sbModel->SetSnapshot(*poi);

    // Create the Background only model form the S+B model
    ModelConfig * bModel = (ModelConfig*) sbModel->Clone();
    bModel->SetName("B Model");      
    poi->setVal(0);
    bModel->SetSnapshot( *poi  );

    //-------------------------------------------------------------

    // Create the AsymptoticCalculator from data,alt model, null model (hypothesis tests using asymptotic properties of likelihood function)
    AsymptoticCalculator  ac(*data, *sbModel, *bModel);
    ac.SetOneSidedDiscovery(true);  // for one-side discovery test
    //ac.SetPrintLevel(-1);  // to suppress print level 

    // Run the calculator
    HypoTestResult * asResult = ac.GetHypoTest();
    asResult->Print();

    // HypoTestInverter
    HypoTestInverter acinverter(ac);

    // Statistical configuration of hypothesis test inverter
    acinverter.SetConfidenceLevel(0.683);
    acinverter.UseCLs(true);

    // Technical configuration of hypothesis test inverter
    acinverter.SetVerbose(false);
    acinverter.SetFixedScan(50,0.0,50.0); // set number of points , xmin and xmax

    // Calculation of limit
    HypoTestInverterResult* acinvresult =  acinverter.GetInterval();

    // Print observed limit
    std::cout << 100*acinverter.ConfidenceLevel() << "%  upper limit : " << acinvresult->UpperLimit() << std::endl;

    //Compute expected limit
    std::cout << "Expected upper limits, using the S+B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << acinvresult->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << acinvresult->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << acinvresult->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << acinvresult->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << acinvresult->GetExpectedUpperLimit(2) << std::endl;

    // Create a CL plot
    TCanvas* acinvcan = new TCanvas();
    HypoTestInverterPlot* acinvplot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",acinvresult);
    acinvplot->Draw("CLb 2CL");  // plot also CLb and CLs+b
    acinvcan->SaveAs("Brazil_plot_asymptotic.png");

    //-------------------------------------------------------------

    std::cout << "\n\nRun now FrequentistCalculator.....\n" << std::endl;
    
    // Create the FrequentistCalculator from data,alt model, null model (frequentist hypothesis test calculators usingtoy data (difference in treatment of nuisanceparameters))

    FrequentistCalculator   fc(*data, *sbModel, *bModel);
    fc.SetToys(500,500);    // 2000 for null (B) and 500 for alt (S+B)

    // Create the test statistics
    ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
    // Use one-sided profile likelihood
    profll.SetOneSidedDiscovery(true);

    // Configure  ToyMCSampler and set the test statistics
    ToyMCSampler *toymcs = (ToyMCSampler*)fc.GetTestStatSampler();
    toymcs->SetTestStatistic(&profll);
    
    if (!sbModel->GetPdf()->canBeExtended())
      toymcs->SetNEventsPerToy(1);
  
    // Run the test
    HypoTestResult * fqResult = fc.GetHypoTest();
    fqResult->Print();

    // Plot test statistic distributions
    TCanvas *can = new TCanvas();
    HypoTestPlot * plot = new HypoTestPlot(*fqResult);
    plot->SetLogYaxis(true);
    plot->Draw();
    can->Draw();  
    can->SaveAs("test_statistic_distributions.png");

    // HypoTestInverter
    HypoTestInverter fcinverter(fc);

    // Statistical configuration of hypothesis test inverter
    fcinverter.SetConfidenceLevel(0.683);
    fcinverter.UseCLs(true);

    // Technical configuration of hypothesis test inverter
    fcinverter.SetVerbose(false);
    fcinverter.SetFixedScan(50,0.0,50.0); // set number of points , xmin and xmax

    // Calculation of limit
    HypoTestInverterResult* fcinvresult =  fcinverter.GetInterval();

    // Print observed limit
    std::cout << 100*fcinverter.ConfidenceLevel() << "%  upper limit : " << fcinvresult->UpperLimit() << std::endl;

    // Compute expected limit
    std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << fcinvresult->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << fcinvresult->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << fcinvresult->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << fcinvresult->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << fcinvresult->GetExpectedUpperLimit(2) << std::endl;

    // Create a CL plot

    TCanvas* fcinvcan = new TCanvas();
    HypoTestInverterPlot* fcinvplot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",fcinvresult);
    fcinvplot->Draw("CLb 2CL");  // plot also CLb and CLs+b
    fcinvcan->Draw();
    fcinvcan->SaveAs("Brazil_plot_frequentist.png");

    //-------------------------------------------------------------
    std::cout << "\n\nRun now ProfileLikelihoodCalculator plc.....\n" << std::endl;

    ProfileLikelihoodCalculator plc(*data, *bModel);
    HypoTestResult *plcResult = plc.GetHypoTest();


    plc.SetConfidenceLevel(0.683);  //68% interval it is equivalent to plc.SetTestSize(0.32);
    //plc.SetConfidenceLevel(0.95); // 95% interval
    LikelihoodInterval * plcinterval = plc.GetInterval();

    const double MLE = poi->getVal();
    const double lowerLimit = plcinterval->LowerLimit(*poi);
    poi->setVal(MLE);
    const double upperLimit = plcinterval->UpperLimit(*poi);
    poi->setVal(MLE);
    std::cout << "68% CL interval: [ " << lowerLimit << " ; " << upperLimit << " ]\n" << std::endl;

    // //One-sided upper-limit
    // LikelihoodInterval * plcinterval = plc.GetInterval();
    // plc.SetTestSize(0.10);
    // const double upperLimit = plcinterval->UpperLimit(*poi);
    // std::cout << "One sided upper limit at 95% CL: "<< upperLimit << std::endl;

    TCanvas *plccan = new TCanvas();
    LikelihoodIntervalPlot plcplot(plcinterval);
    plcplot.SetRange(0,100);
    plcplot.Draw();
    plccan->Draw();
    plccan->SaveAs("Negative_logarithm_of_the_profile _likelihood.png");


    plcResult->Print();
}
