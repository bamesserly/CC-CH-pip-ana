#ifndef PhysicsCalculation_h
#define PhysicsCalculation_h

// Calculate Quantities(Always MeV)
TVector3 AdlerAngle(int RefSystemDef, double dmumom, double dpimom, TVector3 NeuDir, TVector3 MuDir, TVector3 PiDir, double Enu) const;
double CalcQ2(const double Enu, const double Emu, const double Thetamu) const;
double CalcWexp(const double Q2, const double Ehad) const;
double Calcq0(const double Enu, const double Emu) const;
double Calcq3(const double Q2, const double Enu, const double Emu) const;
double Calct(const double epi, const double emu, const double pzpi,
             const double pzmu, const double pxpi, const double pxmu,
             const double pypi, const double pymu) const;

#endif
