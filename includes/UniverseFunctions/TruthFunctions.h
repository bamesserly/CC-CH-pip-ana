#ifndef TruthFunctions_h
#define TruthFunctions_h

virtual double GetEmuTrue() const;
virtual double GetPTmuTrue() const;
virtual double GetPXmuTrue() const;
virtual double GetPYmuTrue() const;
virtual double GetPZmuTrue() const;
virtual double GetPmuTrue() const;
virtual double GetThetamuTrue() const;
virtual double GetThetamuTrueDeg() const;
virtual double GetWexpTrue() const;
virtual double GetWgenie() const;

virtual double GetALRTrue(TruePionIdx) const;
virtual double GetAdlerCosThetaTrue(TruePionIdx) const;
virtual double GetAdlerPhiTrue(TruePionIdx) const;
virtual double GetAllTrackEnergyTrue() const;
virtual double GetPTTrue(TruePionIdx) const;
virtual double GetThetapiTrue(TruePionIdx) const;
virtual double GetThetapiTrueDeg(TruePionIdx) const;
virtual double GetTpiTrue(TruePionIdx) const;
virtual double GetpimuAngleTrue(TruePionIdx) const;
virtual int GetNChargedPionsTrue() const;
virtual int GetPiChargeTrue(TruePionIdx) const;
virtual std::vector<double> GetTpiTrueVec() const;

virtual double GetIntVtxXTrue() const;
virtual double GetIntVtxYTrue() const;
virtual double GetIntVtxZTrue() const;

#endif
