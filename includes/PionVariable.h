#ifndef PionVariable_h
#define PionVariable_h

#include "Variable.h"

template<class T>
class PionVariable : public Variable {
  private:
    typedef std::function<double(const T&)> PointerToTFunction;
    PointerToTFunction pointer_to_GetValue;

  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    PionVariable(const std::string label, const std::string xaxis,
                   const std::string units,
                   const int nbins, const double xmin, const double xmax,
                   PointerToTFunction p = &T::GetDummyVar,
                   const bool is_true = false);

    PionVariable(const std::string label, const std::string xaxis,
                   std::string units, const TArrayD& bins_array,
                   PointerToTFunction p = &T::GetDummyVar,
                   const bool is_true = false);


    //==========================================================================
    // Functions
    //==========================================================================
    // Get the variable's value
    virtual double GetValue (const T& t) const;
};

#endif // PionVariable_h
