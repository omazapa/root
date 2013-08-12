#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

double func(const double* x)
{
	// Rosenbrock
	const double f = 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + 1*(1 - x[0])*(1 - x[0]);
	return f;
}

double gradf(const double* x, unsigned int i)
{
	double df;
	if(i == 0)
		// dRosenbrock/dx[0]
		df = -400*x[0]*(x[1] - x[0]*x[0]) - 2 * (1 - x[0]);
	else
		// dRosenbrock/dx[1]
		df = 200 * (x[1] - x[0]*x[0]);
	return df;
}

int main()
{
	ROOT::Math::GradFunctor gradFunctor(&func, &gradf, 2);
	//ROOT::Math::Functor gradFunctor(&func, 2);

	struct {
		const char* Name;
		const char* Option;
	} minimizerDef[] = {
		{ "Minuit2", "migrad" },
		{ "LVMini", "" },
	};

	for(unsigned int i = 0; i < 2; ++i)
	{
		ROOT::Math::Minimizer* mn =
			 ROOT::Math::Factory::CreateMinimizer(minimizerDef[i].Name, minimizerDef[i].Option);

		std::cout << minimizerDef[i].Name << ":" << std::endl;
		if(!mn)
		{
			std::cout << "  Failed to load " << minimizerDef[i].Name << std::endl;
			continue;
		}

		mn->SetPrintLevel(0);
		mn->SetFunction(gradFunctor);
		mn->SetVariable(0, "x1", -1.2, 0.1);
		mn->SetVariable(1, "x2", 1.0, 0.1);
		mn->Minimize();
		mn->Hesse();

		std::cout << "  Minimization done!" << std::endl;
		std::cout << "  Status=" << mn->Status() << std::endl;

		const double* X = mn->X();
		const double* Err = mn->Errors();
		std::cout << "  Minimum at X=" << X[0] << ", Y=" << X[1] << std::endl;
		std::cout << "  After " << mn->NCalls() << " function evalutions" << std::endl;
		std::cout << "  f=" << mn->MinValue() << std::endl;

		std::cout << "  Errors: " << std::endl;
		std::cout << "    X=" << X[0] << "+/-" << Err[0] << std::endl;
		std::cout << "    Y=" << X[1] << "+/-" << Err[1] << std::endl;

		std::cout << "  Covariance Matrix: " << std::endl;
		std::cout << "    " << mn->CovMatrix(0, 0) << "\t" << mn->CovMatrix(1, 0) << std::endl;
		std::cout << "    " << mn->CovMatrix(0, 1) << "\t" << mn->CovMatrix(1, 1) << std::endl;
	}
}
