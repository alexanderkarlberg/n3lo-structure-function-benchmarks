/// An example in C++ similar to the one found in examples_f90
///
#include <hoppet_v1.h>
#include <apfel/apfelxx.h>
#include<iostream>
#include<cmath>
#include<cstdio>

#include <fstream>

using namespace hoppetv1;

// definition of the initial condition function
void  lha_unpolarized_dummy_pdf(const double & x, const double & Q, double * pdf)
{
  const std::map<int, double> fmap = apfel::QCDEvToPhys(apfel::LHToyPDFs(x, Q));
  for (int i = -6; i <= 6; i++)
    pdf[i+6] = fmap.at(i);
}

//_________________________________________________________________________
std::vector<double> EWCharges(double const& Q)
{
  // Relevant constants
  const double Q2   = Q * Q;
  const double MW2  = apfel::WMass * apfel::WMass;
  const double MZ2  = apfel::ZMass * apfel::ZMass;
  const double s2tw = 1 - MW2 / MZ2;
  const double VD   = - 0.5 + 2 * s2tw / 3;
  const double VU   = + 0.5 - 4 * s2tw / 3;
  const double AD   = - 0.5;
  const double AU   = + 0.5;
  const double Ve   = - 0.5 + 2 * s2tw;
  const double Ae   = - 0.5;
  const std::vector<double> Vq = {VD, VU, VD, VU, VD, VU};
  const std::vector<double> Aq = {AD, AU, AD, AU, AD, AU};
  const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * s2tw * ( 1 - s2tw ) );
  const double PZ2 = PZ * PZ;

  // Build electroweak charges
  std::vector<double> Charges(6, 0.);
  for (auto i = 0; i < 6; i++)
    Charges[i] = apfel::QCh2[i]
      - 2 * apfel::QCh[i] * Vq[i] * Ve * PZ
      + ( Ve * Ve + Ae * Ae ) * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;
  return Charges;
}

int exponent(double const& d) { return (int) floor(log10(fabs(d))); };
std::string expString(double const& d) { return (exponent(d) < 0 ? "-" : "+") + std::to_string((int) fabs(exponent(d))); };
double base(double const& d)  { return d * pow(10, -1.0 * exponent(d)); };

//----------------------------------------------------------------------
int main()
{
  // Evolution parameters
  const double mc = 1.414213563;   
  const double mb = 4.5;
  const double mt = 175.0;
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  const double asQ   = 0.35;
  const double Q0    = sqrt(2.0);
  const double muR_Q = 1.0;

  const int nflav     = -5;
  const int sc_choice = 1;
  const double zmass  = apfel::ZMass;
  const double wmass  = apfel::WMass;
  const double s2tw   = 1 - pow(wmass / zmass, 2);
  const double Ve     = - 0.5 + 2 * s2tw;
  const double Ae     = - 0.5;
  const bool param_coefs = true;

  const double Qmax    = 2.0*13000.0;
  const double Qmin    = 1.0;
  const int order      = -6; 
  const double ymax    = 16.0;
  const double dy      = 0.05;
  const double dlnlnQ  = dy/8.0;

  const bool exact_nfthreshold = true;
  const bool exact_splitting = false;

  const std::vector<int> order_maxv{1, 2, 3, 4};

  for (int k = 0; k < (int) order_maxv.size(); k++)
    {
      const int order_max = order_maxv[k];
      const int nloop     = std::min(order_max, 3);
      std::cout << "# Perturbative order = " << order_max << std::endl;

      // HOPPET
      hoppetSetPoleMassVFN(mc, mb, mt);
      hoppetSetExactDGLAP(exact_nfthreshold, exact_splitting);
      hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme_MSbar);
      hoppetEvolve(asQ, Q0, nloop, muR_Q, lha_unpolarized_dummy_pdf, Q0);

      // APFEL++
      apfel::AlphaQCD a{asQ, Q0, Thresholds, nloop - 1};
      const apfel::TabulateObject<double> Alphas{a, 1000, 0.9, 111, 3};
      const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
      const apfel::Grid g{{apfel::SubGrid{200, 1e-6, 3}, apfel::SubGrid{250, 1e-2, 3}, apfel::SubGrid{200, 6e-1, 3}, apfel::SubGrid{100, 8.5e-1, 3}}};
      const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, Q0, nloop - 1, as);
      const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 500, 1, 110, 3};
      const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x, Q); };
      const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return EWCharges(Q); };
      const auto F2Obj = InitializeF2NCObjectsZM(g, Thresholds);
  
      // Scales
      const std::vector<std::vector<double>> xmuv{{1, 1}, {0.5, 1}, {2, 1}, {1, 0.5}, {1, 2}};

      // Tabulation parameters
      const int nyb = 100;
      //const double ybmin = -3;
      //const double ybmax = 10;
      const double xbmin = 0.00001;
      const double xbmax = 0.95;
      const double ybmin = log(( 1 - xbmax ) / xbmax);
      const double ybmax = log(( 1 - xbmin ) / xbmin);
      const double ybstp = ( ybmax - ybmin ) / nyb;
      const std::vector<double> muv{5, 10, 50, 100};
      double StrFct[14];
  
      std::vector<std::vector<double>> f2h;
      std::vector<std::vector<double>> f2a;
      for (int i = 0; i < (int) xmuv.size(); i++)
	{
	  const double xmuR = xmuv[i][0];
	  const double xmuF = xmuv[i][1];
	  std::cout << "#xmuR = " << xmuR << ", xmuF = " << xmuF << std::endl;

	  // Structure functions
	  hoppetStartStrFctExtended(order_max, nflav, sc_choice, zmass, param_coefs, wmass, zmass);
	  hoppetInitStrFct(order_max, param_coefs, xmuR, xmuF);
	  const auto F2 = BuildStructureFunctions(F2Obj, PDFs, order_max - 1, as, fBq, xmuR, xmuF);

	  std::vector<double> hop;
	  std::vector<double> apf;
	  for (double Q : muv)
	    {
	      const double PZ = pow(Q, 2) / ( pow(Q, 2) + pow(zmass, 2) ) / ( 4 * s2tw * ( 1 - s2tw ) );
	      for (double y = ybmin; y <= 1.0000001 * ybmax; y += ybstp)
		{
		  const double x = 1 / ( 1 + exp(y) );
		  hoppetStrFct(x, Q, xmuR * Q, xmuF * Q, StrFct);
		  hop.push_back(StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * pow(PZ, 2) - StrFct[iF2gZ] * Ve * PZ);
		  apf.push_back(F2.at(0).Evaluate(Q).Evaluate(x));
		}
	    }
	  f2h.push_back(hop);
	  f2a.push_back(apf);
	}

      std::ofstream fout("../plots/F2NC_Scale_Variations_N" + std::to_string(order_max - 1) + "LO.dat");
      fout << std::scientific;
  
      fout << "#x\t\tQ [GeV]\t\tF2[HOPPET]\tF2[APFEL]" << std::endl;
      int i = 0;
      for (double Q : muv)
	{
	  for (double y = ybmin; y <= 1.0000001 * ybmax; y += ybstp)
	    {
	      const double x = 1 / ( 1 + exp(y) );
	      fout << x << "\t" << Q << "\t";
	      for (int j = 0; j < (int) f2h.size(); j++)
		fout << f2h[j][i] << "\t" << f2a[j][i] << "\t";
	      fout << std::endl;
	      i++;
	    }
	  fout << "\n";
	}
      fout.close();
    }
}
