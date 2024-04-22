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
double base(double const& d)  { return d * pow(10, -1.0 * exponent(d)); };

//----------------------------------------------------------------------
int main()
{
  // Evolution parameters
  const double mc = 1.414213563;   
  const double mb = 4.5;
  const double mt = 175.0;
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  const int nflav     = -5;
  const int order_max = 4;
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
  apfel::AlphaQCD a{asQ, Q0, Thresholds, nloop - 1};
  const apfel::TabulateObject<double> Alphas{a, 1000, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
  const apfel::Grid g{{apfel::SubGrid{200, 1e-6, 3}, apfel::SubGrid{250, 1e-2, 3}, apfel::SubGrid{200, 6e-1, 3}, apfel::SubGrid{100, 8.5e-1, 3}}};
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, Q0, nloop - 1, as);
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 500, 1, 1000, 3};
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x, Q); };
  const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return EWCharges(Q); };
  const std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return PVEWCharges(Q); };
  const std::function<std::vector<double>(double const&)> fCKM = [=] (double const&) -> std::vector<double> { return std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1}; };

  // Neutral current structure functions
  const auto F2 = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fBq, xmuR, xmuF);
  const auto FL = BuildStructureFunctions(InitializeFLNCObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fBq, xmuR, xmuF);
  const auto F3 = BuildStructureFunctions(InitializeF3NCObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fDq, xmuR, xmuF);

  // Tabulate
  const apfel::TabulateObject<apfel::Distribution> F2total{[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotal{[&] (double const& Q) -> apfel::Distribution{ return FL.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3total{[&] (double const& Q) -> apfel::Distribution{ return F3.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};

  // Charged current structure functions
  const auto F2p = BuildStructureFunctions(InitializeF2CCPlusObjectsZM(g, Thresholds),  PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
  const auto F2m = BuildStructureFunctions(InitializeF2CCMinusObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
  const auto FLp = BuildStructureFunctions(InitializeFLCCPlusObjectsZM(g, Thresholds),  PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
  const auto FLm = BuildStructureFunctions(InitializeFLCCMinusObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
  const auto F3p = BuildStructureFunctions(InitializeF3CCPlusObjectsZM(g, Thresholds),  PDFs, order_max - 1, as, fCKM, xmuR, xmuF);
  const auto F3m = BuildStructureFunctions(InitializeF3CCMinusObjectsZM(g, Thresholds), PDFs, order_max - 1, as, fCKM, xmuR, xmuF);

  // Tabulate
  const apfel::TabulateObject<apfel::Distribution> F2totalp{[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q) + F2m.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2totalm{[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q) - F2m.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotalp{[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q) + FLm.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotalm{[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q) - FLm.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3totalp{[&] (double const& Q) -> apfel::Distribution { return F3m.at(0).Evaluate(Q) + F3p.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3totalm{[&] (double const& Q) -> apfel::Distribution { return F3m.at(0).Evaluate(Q) - F3p.at(0).Evaluate(Q); }, 500, 1, 400, 3, Thresholds};

  // output the results
  const std::vector<double> xvals{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
  double pdf[13];
  double Q = apfel::ZMass;
  double StrFct[14];
  printf("\n                                Evaluating PDFs and structure functions at Q = %8.4f GeV (HOPPET)\n",Q);
  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon       F1NC        F2NC        F3NC        F1Wp        F1Wm        F2Wp        F2Wm        F3Wp        F3Wm\n");
  for (double x : xvals)
    {
      // propagator
      const double PZ  = pow(Q, 2) / ( pow(Q, 2) + pow(zmass, 2) ) / ( 4 * s2tw * ( 1 - s2tw ) );
      const double PZ2 = PZ * PZ;
      hoppetEval(x, Q, pdf);
      hoppetStrFct(x, Q, xmuR * Q, xmuF * Q,StrFct);
      const std::map<int, double> DistMap = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, Q));
      printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",x,
	     pdf[6+2] - pdf[6-2],
	     pdf[6+1] - pdf[6-1],
	     2 * ( pdf[6-1] + pdf[6-2] ),
	     pdf[6-4] + pdf[6+4],
	     pdf[6+0],
	     StrFct[iF1EM] + StrFct[iF1Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF1gZ] * Ve * PZ,
	     StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF2gZ] * Ve * PZ,
	     2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ,
	     StrFct[iF1Wp],
	     StrFct[iF1Wm],
	     StrFct[iF2Wp],
	     StrFct[iF2Wm],
	     StrFct[iF3Wp],
	     StrFct[iF3Wm]);
    }

  printf("\n                                Evaluating PDFs and structure functions at Q = %8.4f GeV (APFEL++)\n",Q);
  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon       F1NC        F2NC        F3NC        F1Wp        F1Wm        F2Wp        F2Wm        F3Wp        F3Wm\n");
  for (double x : xvals)
    {
      const std::map<int, double> DistMap = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, Q));
      printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",x,
	     DistMap.at(2) - DistMap.at(-2),
	     DistMap.at(1) - DistMap.at(-1),
	     2 * ( DistMap.at(-2) + DistMap.at(-1) ),
	     DistMap.at(4) + DistMap.at(-4),
	     DistMap.at(0),
	     ( F2total.EvaluatexQ(x, Q) - FLtotal.EvaluatexQ(x, Q) ) / 2 / x,
	     F2total.EvaluatexQ(x, Q),
	     F3total.EvaluatexQ(x, Q) / x,
	     ( F2totalp.EvaluatexQ(x, Q) - FLtotalp.EvaluatexQ(x, Q) ) / 2 / x,
	     ( F2totalm.EvaluatexQ(x, Q) - FLtotalm.EvaluatexQ(x, Q) ) / 2 / x,
	     F2totalp.EvaluatexQ(x, Q),
	     F2totalm.EvaluatexQ(x, Q),
	     F3totalp.EvaluatexQ(x, Q) / x,
	     F3totalm.EvaluatexQ(x, Q) / x);
    }
}
