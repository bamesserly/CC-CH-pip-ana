#ifndef HadronVariable2D_cxx
#define HadronVariable2D_cxx

#include "HadronVariable2D.h"


// CTOR -- uniform binning
HadronVariable2D::HadronVariable2D(const std::string labelX, const std::string labelY,
                                   const std::string xaxisX, const std::string xaxisY,
                                   const std::string unitsX, const std::string unitsY,
                                   const int nbinsX, const double xminX, const double xmaxX,
                                   const int nbinsY, const double xminY, const double xmaxY,
                                   PointerToCVUniverseHadronFunction px,
                                   PointerToCVUniverseHadronFunction py,
                                   const bool is_true)
  : Variable2D(labelX, labelY, xaxisX, xaxisY, unitsX, unitsY, nbinsX, xminX, xmaxX, nbinsY, xminY, xmaxY, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), is_true),
    pointer_to_GetHadValueX(px),
    pointer_to_GetHadValueY(py),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
//    m_type(1)
{}

// CTOR -- variable binning
HadronVariable2D::HadronVariable2D(const std::string labelX, const std::string labelY,
                                   const std::string xaxisX, const std::string xaxisY,
                                   const std::string unitsX, const std::string unitsY,
                                   const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                                   PointerToCVUniverseHadronFunction px,
                                   PointerToCVUniverseHadronFunction py,
                                   const bool is_true)
  : Variable2D(labelX, labelY, xaxisX, xaxisY, unitsX, unitsY, bins_arrayX, bins_arrayY, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), is_true),
    pointer_to_GetHadValueX(px),
    pointer_to_GetHadValueY(py),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
  //  m_type(2)
{}

HadronVariable2D::HadronVariable2D(const std::string labelX, const std::string labelY,
                                   const std::string xaxisX, const std::string xaxisY,
                                   const std::string unitsX, const std::string unitsY,
                                   const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                                   PointerToCVUniverseFunction px,
                                   PointerToCVUniverseHadronFunction py,
                                   const bool is_true)
  : Variable2D(labelX, labelY, xaxisX, xaxisY, unitsX, unitsY, bins_arrayX, bins_arrayY, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), is_true),
    m_pointer_to_GetValueX(px),
    pointer_to_GetHadValueY(py),
    pointer_to_GetHadValueX(&CVUniverse::GetDummyHadVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
//    m_type(3)
{}

HadronVariable2D::HadronVariable2D(const std::string labelX, const std::string labelY,
                                   const std::string xaxisX, const std::string xaxisY,
                                   const std::string unitsX, const std::string unitsY,
                                   const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                                   PointerToCVUniverseHadronFunction px,
                                   PointerToCVUniverseFunction py,
                                   const bool is_true, int type )
  : Variable2D(labelX, labelY, xaxisX, xaxisY, unitsX, unitsY, bins_arrayX, bins_arrayY, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), is_true),
    pointer_to_GetHadValueX(px),
    m_pointer_to_GetValueY(py),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar),
    pointer_to_GetHadValueY(&CVUniverse::GetDummyHadVar)
//    m_type(4)
{}
  
// CTOR -- using HadronVariable

HadronVariable2D::HadronVariable2D(const HadronVariable* x,
                                   const HadronVariable* y)
  : Variable2D(x->m_label, y->m_label, x->m_xlabel, y->m_xlabel, x->m_units, y->m_units, x->m_hists.m_bins_array, y->m_hists.m_bins_array, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), x->m_is_true),
    pointer_to_GetHadValueX(x->m_aux_pointer_to_GetHadValue),
    pointer_to_GetHadValueY(y->m_aux_pointer_to_GetHadValue),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
//    m_type(5)
{}

HadronVariable2D::HadronVariable2D(const std::string name,
                                   const HadronVariable* x,
                                   const HadronVariable* y)
  : Variable2D(name, y->m_label, x->m_xlabel, y->m_xlabel, x->m_units, y->m_units, x->m_hists.m_bins_array, y->m_hists.m_bins_array, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), x->m_is_true),
    pointer_to_GetHadValueX(x->m_aux_pointer_to_GetHadValue),
    pointer_to_GetHadValueY(y->m_aux_pointer_to_GetHadValue),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
//    m_type(6)
{}

HadronVariable2D::HadronVariable2D(const HadronVariable* x,
                                   const Variable* y)
  : Variable2D(x->m_label, y->m_label, x->m_xlabel, y->m_xlabel, x->m_units, y->m_units, x->m_hists.m_bins_array, y->m_hists.m_bins_array, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), x->m_is_true),
    pointer_to_GetHadValueX(x->m_aux_pointer_to_GetHadValue),
    m_pointer_to_GetValueY(y->m_aux_pointer_to_GetValue),
    pointer_to_GetHadValueY(&CVUniverse::GetDummyHadVar),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar)
//    m_type(7)
{}

HadronVariable2D::HadronVariable2D(const std::string name,
                                   const HadronVariable* x,
                                   const Variable* y)
  : Variable2D(name, y->m_label, x->m_xlabel, y->m_xlabel, x->m_units, y->m_units, x->m_hists.m_bins_array, y->m_hists.m_bins_array, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), x->m_is_true),
    pointer_to_GetHadValueX(x->m_aux_pointer_to_GetHadValue),
    m_pointer_to_GetValueY(y->m_aux_pointer_to_GetValue),
    pointer_to_GetHadValueY(&CVUniverse::GetDummyHadVar),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar)
//    m_type(8)
{}

HadronVariable2D::HadronVariable2D(const Variable* x,
                                   const HadronVariable* y)
  : Variable2D(x->m_label, y->m_label, x->m_xlabel, y->m_xlabel, x->m_units, y->m_units, x->m_hists.m_bins_array, y->m_hists.m_bins_array, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), x->m_is_true),
    m_pointer_to_GetValueX(x->m_aux_pointer_to_GetValue),
    pointer_to_GetHadValueY(y->m_aux_pointer_to_GetHadValue),
    pointer_to_GetHadValueX(&CVUniverse::GetDummyHadVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
  //  m_type(9)
{}

HadronVariable2D::HadronVariable2D(const std::string name,
                                   const Variable* x,
                                   const HadronVariable* y)
  : Variable2D(name, y->m_label, x->m_xlabel, y->m_xlabel, x->m_units, y->m_units, x->m_hists.m_bins_array, y->m_hists.m_bins_array, PointerToCVUniverseFunction(), PointerToCVUniverseFunction(), x->m_is_true),
    m_pointer_to_GetValueX(x->m_aux_pointer_to_GetValue),
    pointer_to_GetHadValueY(y->m_aux_pointer_to_GetHadValue),
    pointer_to_GetHadValueX(&CVUniverse::GetDummyHadVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar)
//    m_type(10)
{}

// GetValue defines this variable
double HadronVariable2D::GetValueX (const CVUniverse& universe, const int hadron_index) const { 
  if (pointer_to_GetHadValueX(universe, hadron_index) != -9991) return pointer_to_GetHadValueX(universe, hadron_index);
  else return m_pointer_to_GetValueX(universe); 
}
double HadronVariable2D::GetValueY (const CVUniverse& universe, const int hadron_index) const {
  if (pointer_to_GetHadValueY(universe, hadron_index) != -9991) return pointer_to_GetHadValueY(universe, hadron_index);
  else return m_pointer_to_GetValueY(universe);
}

#endif // HadronVariable_h
