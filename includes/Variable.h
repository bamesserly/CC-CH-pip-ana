#ifndef Variable_h
#define Variable_h

#include <functional>
#include "VariableBase.h"

namespace none
{
  struct empty {};
}

#ifndef __CINT__
template<typename T = none::empty>
class Variable : public VariableBase {
 private:
  typedef std::function<double(const CVUniverse&, T)> PointerToCVUniverseFunctionArgs;
  PointerToCVUniverseFunctionArgs m_pointer_to_GetValueArgs;

 public: 
  Variable() {}

  template<typename ...ARGS>
  Variable(PointerToCVUniverseFunctionArgs p, ARGS... args);

  T m_args;

  void SetGetValueArgs(const T& args);

  virtual double GetValue (const CVUniverse& universe) const;
};

#endif // __CINT__

#include "Variable.cxx"

#endif // Variable_h
