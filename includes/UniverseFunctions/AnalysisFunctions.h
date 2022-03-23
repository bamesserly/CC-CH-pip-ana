#ifndef AnalysisFunctions_h
#define AnalysisFunctions_h
// muon
virtual double GetPTmu() const;
virtual double GetThetamuDeg() const;
virtual double GetPXmu() const;
virtual double GetPYmu() const;
virtual double GetPZmu() const;

// event-wide
virtual double GetEhad() const;
virtual double GetEnu() const;
virtual double GetQ2() const;
virtual double Getq0() const;
virtual double Getq3() const;
virtual double GetWexp() const;
virtual int GetNhadrons() const;

// pion
virtual double GetTpi(RecoPionIdx) const;
virtual double GetTpiMBR(RecoPionIdx) const;
virtual double GetPZpi(RecoPionIdx) const;
virtual double GetPXpi(RecoPionIdx) const;
virtual double GetPYpi(RecoPionIdx) const;
virtual double GetPpi(RecoPionIdx) const;
virtual double Gett(RecoPionIdx) const;
virtual double GetThetapi(RecoPionIdx) const;
virtual double GetThetapiDeg(RecoPionIdx) const;
virtual double GetEpi(RecoPionIdx) const;
virtual double GetALR(RecoPionIdx) const;
virtual double GetAdlerCosTheta(RecoPionIdx) const;
virtual double GetAdlerPhi(RecoPionIdx) const;
virtual double GetpimuAngle(RecoPionIdx) const;
virtual double GetPT(RecoPionIdx) const;
#endif
