#ifndef EhadFunctions_h
#define EhadFunctions_h

// Ehad new variables
virtual double GetCalRecoilEnergy() const;
virtual double GetTrackRecoilEnergy() const;
virtual double GetNonCalRecoilEnergy() const;
virtual double GetCalRecoilEnergyNoPi_DefaultSpline() const;
virtual double GetCalRecoilEnergyNoPi_Corrected(const double ecal_nopi) const;
virtual double GetCalRecoilEnergy_DefaultSpline() const;
virtual double GetCalRecoilEnergy_CCPiSpline() const;
virtual double GetCalEpi(RecoPionIdx) const;

// Ehad old variables
virtual double GetCalRecoilEnergy_CCIncSpline() const;
virtual double GetCalRecoilEnergyNoPi_CCIncSpline() const;

// Ehad truth variables
virtual double GetEhadTrue() const;
virtual double GetTpiTrueMatched(RecoPionIdx) const;
virtual double GetEpiTrueMatched(RecoPionIdx) const;
virtual double GetCalRecoilEnergyNoPiTrue() const;

#endif
