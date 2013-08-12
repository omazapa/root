#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFitResult.h>
#include <Fit/FitConfig.h>

double gradient(double *x, double *p, double *u)
{
	u[0] = x[0];
	u[1] = 1.0;
	return p[0]*x[0] + p[1];
}

int main()
{
	TRandom* rand = new TRandom3;
	rand->SetSeed();

	TH1D* h = new TH1D("testHisto", "testHisto", 100, 0.0, 100.0);
	for(int i = 1; i <= h->GetNbinsX(); ++i)
	{
		h->SetBinContent(i, rand->Poisson(i));
		h->SetBinError(i, TMath::Sqrt(h->GetBinContent(i)));
	}

	struct {
		const char* Name;
		const char* Option;
	} minimizerDef[] = {
		{ "Minuit2", "migrad" },
		{ "LVMini", "" },
	};

	for(int i = 0; i < 2; ++i)
	{
		ROOT::Fit::FitConfig::SetDefaultMinimizer(minimizerDef[i].Name, minimizerDef[i].Option);

		TF1* f = new TF1("testFunc", "[0]*x + [1]", 0, 100);
		f->SetGradientFunction(gradient);
		f->SetParName(0, "slope"); f->SetParameter(0, 5.);
		f->SetParName(1, "constant"); f->SetParameter(1, 5.);
		f->CheckGradientFunction();
		TFitResultPtr result = h->Fit(f, "SQG");

		std::cout << minimizerDef[i].Name << " Straight line fit result: " << std::endl;
		std::cout << "  NCalls=" << result->NCalls() << std::endl;
		for(unsigned int i = 0; i < 2; ++i)
			std::cout << "  " << result->ParName(i) << ": " << result->Parameter(i) << " +/- " << result->ParError(i) << std::endl;
		std::cout << "  Correlation=" << 100.*result->GetCorrelationMatrix()(0,1) << "%" << std::endl;
		for(unsigned int i = 0; i < 2; ++i)
			std::cout << "  " << result->ParName(i) << ": " << "GlobalCC=" << 100.*result->GlobalCC(i) << "%" << std::endl;
		//result->PrintOutput();
	}

	return 0;
}
