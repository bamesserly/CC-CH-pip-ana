#ifndef COHDIFFRACTIVESYSTEMATICS_H
#define COHDIFFRACTIVESYSTEMATICS_H

#include "PlotUtils/ChainWrapper.h"
#include "CVUniverse.h"

//We're reweighting coherent events (in CVUniverse) to simulate diffractive events
//Putting a 50% uncertainty on that model choice here
static double fracDiffractiveUnc = 0.5;

//Looking at Alex's coherent xsec in plastic for pion T<500 MeV
// The data w/ error is close to MC 
// The precision on the data is about 20%
// (yes, plastic, but we're assuming its roughly the same for targets)
static double fracCoherentPiUncTracker = 0.2;
static double fracCoherentPiUncNuclear = 0.4;

std::vector<CVUniverse*> GetCohDiffractiveSystematics( PlotUtils::ChainWrapper* chain);
std::map< std::string, std::vector<CVUniverse*> > GetCohDiffractiveSystematicsMap(PlotUtils::ChainWrapper* chain); 
 
//template<class T>
//class DiffractiveUniverse: public T
//XXX Aaron Specific
class DiffractiveUniverse : public CVUniverse
{
  public:
    //DiffractiveUniverse(typename T::config_t chw, double nsigma, double fracDiffUnc = fracDiffractiveUnc );
    DiffractiveUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracDiffUnc = fracDiffractiveUnc );
    
    virtual double GetDiffractiveWeight() const /*override*/;
    
    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  private:
    double m_fracDiffUnc;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//template<class T>
//class CoherentPiUniverse: public T
//XXX Aaron Specific
class CoherentPiPlasticUniverse : public CVUniverse
{
  public:
    //CoherentPiPlasticUniverse(typename T::config_t chw, double nsigma, double fracCohPiUnc = fracCoherentPiUnc );
    CoherentPiPlasticUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc );
    
    virtual double GetCoherentPiWeight( double thpi_true /*deg*/, double tpi_true /*GeV*/ ) const /*override*/;
    
    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  private:
    double m_fracCohPiUnc;
};

class CoherentPiCarbonUniverse : public CVUniverse
{
  public:
    //CoherentPiCarbonUniverse(typename T::config_t chw, double nsigma, double fracCohPiUnc = fracCoherentPiUnc );
    CoherentPiCarbonUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc );
    
    virtual double GetCoherentPiWeight( double thpi_true /*deg*/, double tpi_true /*GeV*/ ) const /*override*/;
    
    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  private:
    double m_fracCohPiUnc;
};

class CoherentPiWaterUniverse : public CVUniverse
{
  public:
    //CoherentPiWaterUniverse(typename T::config_t chw, double nsigma, double fracCohPiUnc = fracCoherentPiUnc );
    CoherentPiWaterUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc );
    
    virtual double GetCoherentPiWeight( double thpi_true /*deg*/, double tpi_true /*GeV*/ ) const /*override*/;
    
    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  private:
    double m_fracCohPiUnc;
};

class CoherentPiIronUniverse : public CVUniverse
{
  public:
    //CoherentPiIronUniverse(typename T::config_t chw, double nsigma, double fracCohPiUnc = fracCoherentPiUnc );
    CoherentPiIronUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc );
    
    virtual double GetCoherentPiWeight( double thpi_true /*deg*/, double tpi_true /*GeV*/ ) const /*override*/;
    
    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  private:
    double m_fracCohPiUnc;
};

class CoherentPiLeadUniverse : public CVUniverse
{
  public:
    //CoherentPiLeadUniverse(typename T::config_t chw, double nsigma, double fracCohPiUnc = fracCoherentPiUnc );
    CoherentPiLeadUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc );
    
    virtual double GetCoherentPiWeight( double thpi_true /*deg*/, double tpi_true /*GeV*/ ) const /*override*/;
    
    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  private:
    double m_fracCohPiUnc;
};


#endif // COHDIFFRACTIVESYSTEMATICS_H

