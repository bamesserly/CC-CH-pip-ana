#ifndef Variable_cxx
#define Variable_cxx

#include "Variable.h"

#ifndef __CINT__ // CINT doesn't know about std::function

template<typename T>
template<typename ...ARGS>
Variable<T>::Variable(PointerToCVUniverseFunctionArgs p, ARGS... args) 
  : VariableBase(args...), m_pointer_to_GetValueArgs(p), m_args()
{}

template<typename T>
void Variable<T>::SetGetValueArgs(const T& args) { m_args = args; }

template<typename T>
double Variable<T>::GetValue (const CVUniverse& universe) const {
  return m_pointer_to_GetValueArgs(universe, m_args);
}

#endif // __CINT__

#endif // Variable_cxx
