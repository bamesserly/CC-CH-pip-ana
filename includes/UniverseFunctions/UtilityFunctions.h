#ifndef UtilityFunctions_h
#define UtilityFunctions_h

// Variables used for cuts
// As well as retired variables from past studies
double GetTpiFResidual(const int hadron, const bool MBR = false) const;
double GetWexpFResidual() const;
virtual bool IsInHexagon(double x, double y, double apothem) const;
virtual bool IsInPlastic() const;
virtual bool leftlinesCut(const double a, const double x, const double y) const;
virtual bool rightlinesCut(const double a, const double x,
                           const double y) const;
virtual double GetEmichel(RecoPionIdx) const;
virtual double GetEnode0(RecoPionIdx) const;
virtual double GetEnode01(RecoPionIdx) const;
virtual double GetEnode1(RecoPionIdx) const;
virtual double GetEnode2(RecoPionIdx) const;
virtual double GetEnode3(RecoPionIdx) const;
virtual double GetEnode4(RecoPionIdx) const;
virtual double GetEnode5(RecoPionIdx) const;
virtual double GetFitVtxX() const;
virtual double GetFitVtxY() const;
virtual double GetFitVtxZ() const;
virtual double GetLLRScore(RecoPionIdx) const;
virtual double GetLargestIsoProngSep() const;
virtual double GetLargestPrimProngSep() const;
virtual double GetdEdxScore(RecoPionIdx) const;
virtual int GetNAnchoredLongTracks() const;
virtual int GetNAnchoredShortTracks() const;
virtual int GetNIsoProngs() const;
virtual int GetNNodes(RecoPionIdx) const;
virtual int GetNhadrons() const;
virtual int GetTrackReconstructionMethod(RecoPionIdx) const;

#endif
