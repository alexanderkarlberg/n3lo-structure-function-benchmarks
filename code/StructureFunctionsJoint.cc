/// An example in C++ similar to the one found in examples_f90
///
#include <hoppet_v1.h>
#include <apfel/apfelxx.h>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace hoppetv1;
// definition of the initial condition function
void  lha_unpolarized_dummy_pdf(const double & x, const double & Q, double * pdf)
{
  for (auto const& f : apfel::QCDEvToPhys(apfel::LHToyPDFs(x, Q)))
    pdf[f.first + 6] = f.second;
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

//_________________________________________________________________________
std::vector<double> PVEWCharges(double const& Q)
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
    Charges[i] = - 2 * apfel::QCh[i] * Aq[i] * Ae * PZ + 4 * Vq[i] * Aq[i] * Ve * Ae * PZ2;
  return Charges;
}

int exponent(double const& d) { return (int) floor(log10(fabs(d))); };
std::string expString(double const& d) { return (exponent(d) < 0 ? "-" : "+") + std::to_string((int) fabs(exponent(d))); };
double base(double const& d) { return d * pow(10, -1.0 * exponent(d)); };

//----------------------------------------------------------------------
int main()
{
  // Evolution parameters
  const double mc = 1.414213563;
  const double mb = 4.5;
  const double mt = 175.0;
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  const std::vector<int> vorder{1, 2, 3, 4};
  const std::vector<std::string> vorder_str{"LO", "NLO", "NNLO", "N$^{3}$LO"};
  const std::vector<std::string> vorderlog_str{"LL", "NLL", "NNLL", "N$^{3}$LL"};

  for (int order_max : vorder)
    {
      const int nflav     = -5;
      const int sc_choice = 1 ;
      const double zmass  = apfel::ZMass;
      const double wmass  = apfel::WMass;
      const double s2tw   = 1 - pow(wmass / zmass, 2);
      const double Ve     = - 0.5 + 2 * s2tw;
      const double Ae     = - 0.5;
      const bool param_coefs = true;

      const double Qmax    = 13000.0;
      const double Qmin    = 1.0;
      const int order      = -6; 
      const double ymax    = 16.0;
      const double dy      = 0.05;
      const double dlnlnQ  = dy/8.0;
      const int    nloop   = std::min(order_max, 3);
      const double xmuR    = 1.0;
      const double xmuF    = 1.0;
      const double minQval = std::min(xmuF * Qmin, Qmin);
      const double maxQval = std::max(xmuF * Qmax, Qmax);

      const double asQ   = 0.35;
      const double Q0    = sqrt(2.0);
      const double muR_Q = 1.0;

      const bool exact_nfthreshold = true;
      const bool exact_splitting = false;
  
      // HOPPET
      hoppetSetPoleMassVFN(mc,mb,mt);
      hoppetSetExactDGLAP(exact_nfthreshold,exact_splitting);
      hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,order,factscheme_MSbar);
      hoppetStartStrFctExtended(order_max, nflav,sc_choice,zmass,param_coefs,wmass,zmass);
      hoppetEvolve(asQ, Q0, nloop, muR_Q, lha_unpolarized_dummy_pdf, Q0);
      hoppetInitStrFct(order_max,param_coefs,xmuR,xmuF);

      // APFEL++
      apfel::Banner();
      apfel::AlphaQCD a{asQ, Q0, Thresholds, nloop - 1};
      const apfel::TabulateObject<double> Alphas{a, 1000, 0.9, 1001, 3};
      const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
      const apfel::Grid g{{apfel::SubGrid{200, 1e-6, 5}, apfel::SubGrid{250, 1e-2, 5}, apfel::SubGrid{200, 6e-1, 5}, apfel::SubGrid{100, 8.5e-1, 5}}};
      const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, Q0, nloop - 1, as);
      const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 1000, 1, 800, 3};
      const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x, Q); };
      const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return EWCharges(Q); };
      const std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return PVEWCharges(Q); };
      const std::function<std::vector<double>(double const&)> fCKM = [=] (double const&) -> std::vector<double> { return std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1}; };

      // Neutral current structure functions
      const auto F2 = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fBq, xmuR, xmuF);
      const auto FL = BuildStructureFunctions(InitializeFLNCObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fBq, xmuR, xmuF);
      const auto F3 = BuildStructureFunctions(InitializeF3NCObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fDq, xmuR, xmuF);

      // Tabulate
      const apfel::TabulateObject<apfel::Distribution> F2total{[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLtotal{[&] (double const& Q) -> apfel::Distribution{ return FL.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3total{[&] (double const& Q) -> apfel::Distribution{ return F3.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};

      // Charged current structure functions
      const auto F2p = BuildStructureFunctions(InitializeF2CCPlusObjectsZM(g, Thresholds),  PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
      const auto F2m = BuildStructureFunctions(InitializeF2CCMinusObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
      const auto FLp = BuildStructureFunctions(InitializeFLCCPlusObjectsZM(g, Thresholds),  PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
      const auto FLm = BuildStructureFunctions(InitializeFLCCMinusObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
      const auto F3p = BuildStructureFunctions(InitializeF3CCPlusObjectsZM(g, Thresholds),  PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
      const auto F3m = BuildStructureFunctions(InitializeF3CCMinusObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fCKM, xmuR, xmuF);

      // Tabulate
      const apfel::TabulateObject<apfel::Distribution> F2totalp{[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q) + F2m.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2totalm{[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q) - F2m.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLtotalp{[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q) + FLm.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLtotalm{[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q) - FLm.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3totalp{[&] (double const& Q) -> apfel::Distribution { return F3m.at(0).Evaluate(Q) + F3p.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3totalm{[&] (double const& Q) -> apfel::Distribution { return F3m.at(0).Evaluate(Q) - F3p.at(0).Evaluate(Q); }, 500, 1, 600, 3, Thresholds};

      // Number for the plots
      const int nyb = 500;
      //const double ybmin = -3;
      //const double ybmax = 10;
      const double xbmin = 0.00001;
      const double xbmax = 0.95;
      const double ybmin = log(( 1 - xbmax ) / xbmax);
      const double ybmax = log(( 1 - xbmin ) / xbmin);
      const double ybstp = ( ybmax - ybmin ) / nyb;
      const std::vector<double> muv{2, 5, 10, 50, 100};

      std::ofstream fout("../plots/StructureFunctions_N" + std::to_string(order_max - 1) + "LO.dat");
      fout << std::scientific;
      double pdf[13];
      double StrFct[14];
      for (double Q : muv)
	{
	  const double PZ  = pow(Q, 2) / ( pow(Q, 2) + pow(zmass, 2) ) / ( 4 * s2tw * ( 1 - s2tw ) );
	  const double PZ2 = PZ * PZ;
	  fout << "#x\t\t"
	       << "Q [GeV]\t\t"
	       << "F1NC [HOPPET]\t"
	       << "F1NC [APFEL]\t"
	       << "F2NC [HOPPET]\t"
	       << "F2NC [APFEL]\t"
	       << "F3NC [HOPPET]\t"
	       << "F3NC [APFEL]\t"
	       << "F1Wp [HOPPET]\t"
	       << "F1Wp [APFEL]\t"
	       << "F1Wm [HOPPET]\t"
	       << "F1Wm [APFEL]\t"
	       << "F2Wp [HOPPET]\t"
	       << "F2Wp [APFEL]\t"
	       << "F2Wm [HOPPET]\t"
	       << "F2Wm [APFEL]\t"
	       << "F3Wp [HOPPET]\t"
	       << "F3Wp [APFEL]\t"
	       << "F3Wm [HOPPET]\t"
	       << "F3Wm [APFEL]\t"
	       << "XNCp [HOPPET]\t"
	       << "XNCp [APFEL]\t"
	       << "XNCm [HOPPET]\t"
	       << "XNCm [APFEL]\t"
	       << "XCCp [HOPPET]\t"
	       << "XCCp [APFEL]\t"
	       << "XCCm [HOPPET]\t"
	       << "XCCm [APFEL]\t"
	       << "\n";
	  for (double y = ybmin; y <= 1.0000001 * ybmax; y += ybstp)
	    {
	      const double x = 1 / ( 1 + exp(y) );
	      hoppetStrFct(x, Q, xmuR * Q, xmuF * Q, StrFct);

	      // Construct structure functions
	      const double F1NCh = StrFct[iF1EM] + StrFct[iF1Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF1gZ] * Ve * PZ;
	      const double F2NCh = StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF2gZ] * Ve * PZ;
	      const double F3NCh = 2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ;
	      const double FLNCh = F2NCh - 2 * x * F1NCh;

	      const double FLNCa = FLtotal.EvaluatexQ(x, Q);
	      const double F2NCa = F2total.EvaluatexQ(x, Q);
	      const double F3NCa = F3total.EvaluatexQ(x, Q) / x;
	      const double F1NCa = ( F2NCa - FLNCa ) / 2 / x;

	      const double F1CCph = StrFct[iF1Wp];
	      const double F2CCph = StrFct[iF2Wp];
	      const double F3CCph = StrFct[iF3Wp];
	      const double FLCCph = F2CCph - 2 * x * F1CCph;

	      const double FLCCpa = FLtotalp.EvaluatexQ(x, Q);
	      const double F2CCpa = F2totalp.EvaluatexQ(x, Q);
	      const double F3CCpa = F3totalp.EvaluatexQ(x, Q) / x;
	      const double F1CCpa = ( F2CCpa - FLCCpa ) / 2 / x;

	      const double F1CCmh = StrFct[iF1Wm];
	      const double F2CCmh = StrFct[iF2Wm];
	      const double F3CCmh = StrFct[iF3Wm];
	      const double FLCCmh = F2CCmh - 2 * x * F1CCmh;

	      const double FLCCma = FLtotalm.EvaluatexQ(x, Q);
	      const double F2CCma = F2totalm.EvaluatexQ(x, Q);
	      const double F3CCma = F3totalm.EvaluatexQ(x, Q) / x;
	      const double F1CCma = ( F2CCma - FLCCma ) / 2 / x;

	      // Construct reduced cross sections assuming a c.m.e. typical
	      // of HERA
	      const double Vs = 320;
	      const double ye = pow(Q / Vs, 2) / x;
	      const double yp = 1 + pow(1 - ye, 2);
	      const double ym = 1 - pow(1 - ye, 2);

	      const double RedXSecNCph = F2NCh - pow(ye, 2) / yp * FLNCh - ym / yp * x * F3NCh;
	      const double RedXSecNCmh = F2NCh - pow(ye, 2) / yp * FLNCh + ym / yp * x * F3NCh;

	      const double RedXSecNCpa = F2NCa - pow(ye, 2) / yp * FLNCa - ym / yp * x * F3NCa;
	      const double RedXSecNCma = F2NCa - pow(ye, 2) / yp * FLNCa + ym / yp * x * F3NCa;

	      const double RedXSecCCph = yp * F2CCph - pow(ye, 2) * FLCCph - ym * x * F3CCph;
	      const double RedXSecCCmh = yp * F2CCmh - pow(ye, 2) * FLCCmh + ym * x * F3CCmh;

	      const double RedXSecCCpa = yp * F2CCpa - pow(ye, 2) * FLCCpa - ym * x * F3CCpa;
	      const double RedXSecCCma = yp * F2CCma - pow(ye, 2) * FLCCma + ym * x * F3CCma;

	      fout << x << "\t" << Q << "\t"
		   << F1NCh  << "\t" << F1NCa  << "\t"
		   << F2NCh  << "\t" << F2NCa  << "\t"
		   << F3NCh  << "\t" << F3NCa  << "\t"
		   << F1CCph << "\t" << F1CCpa << "\t"
		   << F1CCmh << "\t" << F1CCma << "\t"
		   << F2CCph << "\t" << F2CCpa << "\t"
		   << F2CCmh << "\t" << F2CCma << "\t"
		   << F3CCph << "\t" << F3CCpa << "\t"
		   << F3CCmh << "\t" << F3CCma << "\t"
		   << (ye > 1 ? std::nan(0) : RedXSecNCph) << "\t" << (ye > 1 ? std::nan(0) : RedXSecNCpa) << "\t"
		   << (ye > 1 ? std::nan(0) : RedXSecNCmh) << "\t" << (ye > 1 ? std::nan(0) : RedXSecNCma) << "\t"
		   << (ye > 1 ? std::nan(0) : RedXSecCCph) << "\t" << (ye > 1 ? std::nan(0) : RedXSecCCpa) << "\t"
		   << (ye > 1 ? std::nan(0) : RedXSecCCmh) << "\t" << (ye > 1 ? std::nan(0) : RedXSecCCma)
		   << std::endl;
	    }
	  fout << "\n";
	}
      fout.close();

      // Output the results
      // Create latex tables with both HOPPET and APFEL
      FILE* pFileh = fopen(("../tables/table_N" + std::to_string(order_max - 1) + "LO_HOPPET.tex").c_str(), "w");
      FILE* pFilea = fopen(("../tables/table_N" + std::to_string(order_max - 1) + "LO_APFEL.tex").c_str(), "w");
      std::vector<double> xvals{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
      hoppetInitStrFct(order_max,param_coefs,xmuR,xmuF);
      for (double Q : std::vector<double>{2, 50, 100})
	{
	  const double PZ  = pow(Q, 2) / ( pow(Q, 2) + pow(zmass, 2) ) / ( 4 * s2tw * ( 1 - s2tw ) );
	  const double PZ2 = PZ * PZ;

	  printf("                                Evaluating PDFs and structure functions at Q = %8.4f GeV (|HOPPET / APFEL++ - 1|)\n",Q);
	  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon       F1NC        F2NC        F3NC        F1Wp        F1Wm        F2Wp        F2Wm        F3Wp        F3Wm\n");

	  fprintf(pFileh, "\\begin{table}[h]\n");
	  fprintf(pFileh, "\\begin{adjustbox}{width=\\textwidth}\n");
	  fprintf(pFileh, "\\begin{tabular}{|c||c|c|c|c|c|c|c|c|c|}\n");
	  fprintf(pFileh, "\\hline\n");
	  fprintf(pFileh, "$x_{\\rm B}$ & $F_1^{\\rm NC}$ & $F_2^{\\rm NC}$ & $F_3^{\\rm NC}$ & $F_1^{W^+}$ & $F_2^{W^+}$ & $F_3^{W^+}$ & $F_1^{W^-}$ & $F_2^{W^-}$ & $F_3^{W^-}$ \\\\\n");
	  fprintf(pFileh, "\\hline\n");

	  fprintf(pFilea, "\\begin{table}[h]\n");
	  fprintf(pFilea, "\\begin{adjustbox}{width=\\textwidth}\n");
	  fprintf(pFilea, "\\begin{tabular}{|c||c|c|c|c|c|c|c|c|c|}\n");
	  fprintf(pFilea, "\\hline\n");
	  fprintf(pFilea, "$x_{\\rm B}$ & $F_1^{\\rm NC}$ & $F_2^{\\rm NC}$ & $F_3^{\\rm NC}$ & $F_1^{W^+}$ & $F_2^{W^+}$ & $F_3^{W^+}$ & $F_1^{W^-}$ & $F_2^{W^-}$ & $F_3^{W^-}$ \\\\\n");
	  fprintf(pFilea, "\\hline\n");

	  for (double x : xvals)
	    {
	      hoppetEval(x, Q, pdf);
	      hoppetStrFct(x, Q, xmuR * Q, xmuF * Q, StrFct);
	      const std::map<int, double> DistMap = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, Q));

	      printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",x,
		     std::abs(( pdf[6+2] - pdf[6-2] ) / ( DistMap.at(2) - DistMap.at(-2) ) - 1.),
		     std::abs(( pdf[6+1] - pdf[6-1] ) / ( DistMap.at(1) - DistMap.at(-1) ) - 1.),
		     std::abs(( pdf[6-1] + pdf[6-2] ) / ( DistMap.at(-2) + DistMap.at(-1) ) - 1.),
		     std::abs(( pdf[6-4] + pdf[6+4] ) / ( DistMap.at(4) + DistMap.at(-4) ) - 1.),
		     std::abs(pdf[6+0] / DistMap.at(0) - 1.),
		     std::abs(( StrFct[iF1EM] + StrFct[iF1Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF1gZ] * Ve * PZ ) / ( ( F2total.EvaluatexQ(x, Q) - FLtotal.EvaluatexQ(x, Q) ) / 2 / x ) - 1.),
		     std::abs(( StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF2gZ] * Ve * PZ ) / F2total.EvaluatexQ(x, Q) - 1.),
		     std::abs(( 2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ ) / ( F3total.EvaluatexQ(x, Q) / x ) - 1.),
		     std::abs(StrFct[iF1Wp] / ( ( F2totalp.EvaluatexQ(x, Q) - FLtotalp.EvaluatexQ(x, Q) ) / 2 / x ) - 1.),
		     std::abs(StrFct[iF1Wm] / ( ( F2totalm.EvaluatexQ(x, Q) - FLtotalm.EvaluatexQ(x, Q) ) / 2 / x ) - 1.),
		     std::abs(StrFct[iF2Wp] / F2totalp.EvaluatexQ(x, Q) - 1.),
		     std::abs(StrFct[iF2Wm] / F2totalm.EvaluatexQ(x, Q) - 1.),
		     std::abs(StrFct[iF3Wp] / ( F3totalp.EvaluatexQ(x, Q) / x ) - 1.),
		     std::abs(StrFct[iF3Wm] / ( F3totalm.EvaluatexQ(x, Q) / x ) - 1.));

	      fprintf(pFileh, "$%4.1f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ \\\\\n",
		      base(x), expString(x).c_str(),
		      base(StrFct[iF1EM] + StrFct[iF1Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF1gZ] * Ve * PZ), expString(StrFct[iF1EM] + StrFct[iF1Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF1gZ] * Ve * PZ).c_str(),
		      base(StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF2gZ] * Ve * PZ), expString(StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF2gZ] * Ve * PZ).c_str(),
		      base(2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ), expString(2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ).c_str(),
		      base(StrFct[iF1Wp]), expString(StrFct[iF1Wp]).c_str(),
		      base(StrFct[iF1Wm]), expString(StrFct[iF1Wm]).c_str(),
		      base(StrFct[iF2Wp]), expString(StrFct[iF2Wp]).c_str(),
		      base(StrFct[iF2Wm]), expString(StrFct[iF2Wm]).c_str(),
		      base(StrFct[iF3Wp]), expString(StrFct[iF3Wp]).c_str(),
		      base(StrFct[iF3Wm]), expString(StrFct[iF3Wm]).c_str());

	      fprintf(pFilea, "$%4.1f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ & $%7.4f^{%s}$ \\\\\n",
		      base(x), expString(x).c_str(),
		      base(( F2total.EvaluatexQ(x, Q) - FLtotal.EvaluatexQ(x, Q) ) / 2 / x), expString(( F2total.EvaluatexQ(x, Q) - FLtotal.EvaluatexQ(x, Q) ) / 2 / x).c_str(),
		      base(F2total.EvaluatexQ(x, Q)), expString(F2total.EvaluatexQ(x, Q)).c_str(),
		      base(F3total.EvaluatexQ(x, Q) / x), expString(F3total.EvaluatexQ(x, Q) / x).c_str(),
		      base(( F2totalp.EvaluatexQ(x, Q) - FLtotalp.EvaluatexQ(x, Q) ) / 2 / x), expString(( F2totalp.EvaluatexQ(x, Q) - FLtotalp.EvaluatexQ(x, Q) ) / 2 / x).c_str(),
		      base(( F2totalm.EvaluatexQ(x, Q) - FLtotalm.EvaluatexQ(x, Q) ) / 2 / x), expString(( F2totalm.EvaluatexQ(x, Q) - FLtotalm.EvaluatexQ(x, Q) ) / 2 / x).c_str(),
		      base(F2totalp.EvaluatexQ(x, Q)), expString(F2totalp.EvaluatexQ(x, Q)).c_str(),
		      base(F2totalm.EvaluatexQ(x, Q)), expString(F2totalm.EvaluatexQ(x, Q)).c_str(),
		      base(F3totalp.EvaluatexQ(x, Q) / x), expString(F3totalp.EvaluatexQ(x, Q) / x).c_str(),
		      base(F3totalm.EvaluatexQ(x, Q) / x), expString(F3totalm.EvaluatexQ(x, Q) / x).c_str());
	    }
	  printf("\n");

	  fprintf(pFileh, "\\hline\n");
	  fprintf(pFileh, "\\end{tabular}\n");
	  fprintf(pFileh, "\\end{adjustbox}");
	  fprintf(pFileh, "%s", ("\\caption{" + vorder_str[order_max - 1] + " stucture functions with " + vorderlog_str[nloop - 1] + " evolution at $Q = " + std::to_string((int) Q) + "$ GeV.}\n").c_str());
	  fprintf(pFileh, "%s", ("\\label{tab:N" + std::to_string(order_max - 1) + "LO-Q" + std::to_string((int) Q) + "}\n").c_str());
	  fprintf(pFileh, "\\end{table}\n\n\n");

	  fprintf(pFilea, "\\hline\n");
	  fprintf(pFilea, "\\end{tabular}\n");
	  fprintf(pFilea, "\\end{adjustbox}");
	  fprintf(pFilea, "%s", ("\\caption{" + vorder_str[order_max - 1] + " stucture functions with " + vorderlog_str[nloop - 1] + " evolution at $Q = " + std::to_string((int) Q) + "$ GeV.}\n").c_str());
	  fprintf(pFilea, "%s", ("\\label{tab:N" + std::to_string(order_max - 1) + "LO-Q" + std::to_string((int) Q) + "}\n").c_str());
	  fprintf(pFilea, "\\end{table}\n\n\n");
	}
      fclose(pFileh);
      fclose(pFilea);
    }
}
