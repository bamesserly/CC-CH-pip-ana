#ifndef HadronVariable_h
#define HadronVariable_h

#include "Variable.h"

class HadronVariable : public Variable {
  private:
    typedef std::function<double(const CVUniverse&, int)> PointerToCVUniverseHadronFunction ;
    PointerToCVUniverseHadronFunction pointer_to_GetHadValue;

  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    HadronVariable(const std::string label, const std::string xaxis,
                   const std::string units,
                   const int nbins, const double xmin, const double xmax,
                   PointerToCVUniverseHadronFunction p = &CVUniverse::GetDummyHadVar,
                   const bool is_true = false);

    HadronVariable(const std::string label, const std::string xaxis,
                   std::string units, const TArrayD& bins_array,
                   PointerToCVUniverseHadronFunction p = &CVUniverse::GetDummyHadVar,
                   const bool is_true = false);


    //==========================================================================
    // Functions
    //==========================================================================
    // Get the variable's value
    virtual double GetValue (const CVUniverse& universe, const int hadron_index) const;
};

#endif // HadronVariable_h
