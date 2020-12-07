#ifndef WSidebandFitter_h
#define WSidebandFitter_h

#include "PlotUtils/MnvH1D.h"
#include "CVUniverse.h"
#include "Histograms.h"
#include "TFractionFitter.h"

// Alternative fitting method using TFractionFitter
void FractionFit(TH1D* data_hist,              TObjArray* mc_templates,
                 std::vector<double>& weights, std::vector<double>& errs,
                 bool be_verbose,
                 double& midscale, double& hiscale);

static const int kNCoefficients = 3; // For chi2 method
static const int kLoWParamId    = 0;
static const int kMidWParamId   = 1;
static const int kHiWParamId    = 2;
static const int kSigParamId    = 3; // For FractionFitter method


class WSidebandFitter {
  public:
    WSidebandFitter();
    WSidebandFitter(const CVUniverse& universe, Histograms hists, const double pot_scale);
    ~WSidebandFitter();

    static TH1D* m_data;
    static TH1D* m_sig;
    static TH1D* m_loW;
    static TH1D* m_midW;
    static TH1D* m_hiW;

    static double MinimizerFunc(const double* par); // Function that gets passed to minimizer, must have this signature.
    void Fit(); // Do the thing

    static double m_pot_scale;
    static double m_chi2;  // determined in fit
    static int m_ndf;      // determined in fit
    double m_fit_scale[kNCoefficients];     // results of the fit. 0th index is loW, 1st index is hiW.
    double m_fit_scale_err[kNCoefficients]; // errors on results of the fit
    double m_fit_min;                       // returned from the fit

  private:
};

#endif // WSidebandFitter_h
