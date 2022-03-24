#ifndef EhadFunctions_h
#define EhadFunctions_h

// Ehad new variables
virtual double GetCalEpi(RecoPionIdx) const;
virtual double GetCalRecoilEnergy() const;
virtual double GetCalRecoilEnergyNoPi_Corrected(const double ecal_nopi) const;
virtual double GetCalRecoilEnergyNoPi_DefaultSpline() const;
virtual double GetCalRecoilEnergy_CCPiSpline() const;
virtual double GetCalRecoilEnergy_DefaultSpline() const;
virtual double GetNonCalRecoilEnergy() const;
virtual double GetTrackRecoilEnergy() const;

// Ehad old variables
virtual double GetCalRecoilEnergyNoPi_CCIncSpline() const;
virtual double GetCalRecoilEnergy_CCIncSpline() const;

// Ehad truth variables
virtual double GetCalRecoilEnergyNoPiTrue() const;
virtual double GetEhadTrue() const;
virtual double GetEpiTrueMatched(RecoPionIdx) const;
virtual double GetTpiTrueMatched(RecoPionIdx) const;

#endif
