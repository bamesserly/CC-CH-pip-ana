#ifndef HadronVariable_cxx
#define HadronVariable_cxx

#include "HadronVariable.h"


// CTOR -- uniform binning
HadronVariable::HadronVariable(const std::string label, const std::string xaxis,
                               const std::string units,
                               const int nbins, const double xmin, const double xmax,
                               PointerToCVUniverseHadronFunction p,
                               const bool is_true)
  : Variable(label, xaxis, units, nbins, xmin, xmax, PointerToCVUniverseFunction(), is_true),
    pointer_to_GetHadValue(p),
    m_aux_pointer_to_GetHadValue(pointer_to_GetHadValue)
{}


// CTOR -- variable binning
HadronVariable::HadronVariable(const std::string label, const std::string xaxis,
                               const std::string units, const TArrayD& bins_array,
                               PointerToCVUniverseHadronFunction p,
                               const bool is_true)
  : Variable(label, xaxis, units, bins_array, PointerToCVUniverseFunction(), is_true),
    pointer_to_GetHadValue(p),
    m_aux_pointer_to_GetHadValue(pointer_to_GetHadValue)
{}
  

// GetValue defines this variable
double HadronVariable::GetValue (const CVUniverse& universe, const int hadron_index) const { 
  return pointer_to_GetHadValue(universe, hadron_index); 
}


#endif // HadronVariable_h
