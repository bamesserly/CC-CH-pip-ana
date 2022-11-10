#ifndef HadronVariable2D_h
#define HadronVariable2D_h

#include "Variable2D.h"
#include "HadronVariable.h"
class HadronVariable;
class Variable;

class HadronVariable2D : public Variable2D {
  private:
    typedef std::function<double(const CVUniverse&, int)> PointerToCVUniverseHadronFunction ;
    PointerToCVUniverseHadronFunction pointer_to_GetHadValueX;
    PointerToCVUniverseHadronFunction pointer_to_GetHadValueY;
    PointerToCVUniverseFunction m_pointer_to_GetValueX;
    PointerToCVUniverseFunction m_pointer_to_GetValueY;

  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    HadronVariable2D(const std::string labelX, const std::string labelY,
                     const std::string xaxisX, const std::string xaxisY,
                     const std::string unitsX, const std::string unitsY,
                     const int nbinsX, const double xminX, const double xmaxX,
                     const int nbinsY, const double xminY, const double xmaxY,
                     PointerToCVUniverseHadronFunction px = &CVUniverse::GetDummyHadVar,
                     PointerToCVUniverseHadronFunction py = &CVUniverse::GetDummyHadVar,
                     const bool is_true = false);

    HadronVariable2D(const std::string labelX, const std::string labelY,
                     const std::string xaxisX, const std::string xaxisY,
                     const std::string unitsX, const std::string unitsY,
                     const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                     PointerToCVUniverseHadronFunction px = &CVUniverse::GetDummyHadVar,
                     PointerToCVUniverseHadronFunction py = &CVUniverse::GetDummyHadVar,
                     const bool is_true = false);

    HadronVariable2D(const std::string labelX, const std::string labelY,
                     const std::string xaxisX, const std::string xaxisY,
                     const std::string unitsX, const std::string unitsY,
                     const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                     PointerToCVUniverseFunction px = &CVUniverse::GetDummyVar,
                     PointerToCVUniverseHadronFunction py = &CVUniverse::GetDummyHadVar,
                     const bool is_true = false);

    HadronVariable2D(const std::string labelX, const std::string labelY,
                     const std::string xaxisX, const std::string xaxisY,
                     const std::string unitsX, const std::string unitsY,
                     const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                     PointerToCVUniverseHadronFunction px = &CVUniverse::GetDummyHadVar,
                     PointerToCVUniverseFunction py = &CVUniverse::GetDummyVar,
                     const bool is_true = false, int type = -1);

    HadronVariable2D(const HadronVariable*,
                     const HadronVariable*);

    HadronVariable2D(const std::string name,
                     const HadronVariable*,
                     const HadronVariable*);

    HadronVariable2D(const HadronVariable*,
                     const Variable*);

    HadronVariable2D(const std::string name,
                     const HadronVariable*,
                     const Variable*);

    HadronVariable2D(const Variable*,
                     const HadronVariable*);

    HadronVariable2D(const std::string name,
                     const Variable*,
                     const HadronVariable*);


//    int m_type;

    //==========================================================================
    // Functions
    //==========================================================================
    // Get the variable's value
    virtual double GetValueX (const CVUniverse& universe, const int hadron_index) const;

    virtual double GetValueY (const CVUniverse& universe, const int hadron_index) const;
};

#endif // HadronVariable2D_h
