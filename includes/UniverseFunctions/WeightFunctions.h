#ifndef WeightFunctions_h
#define WeightFunctions_h

// Get Weight
virtual double GetDiffractiveWeight() const;
virtual double GetWeight() const;
// Warping
virtual double GetAnisoDeltaDecayWarpWeight() const;
virtual double GetGenieWarpWeight() const;
virtual double GetLowQ2PiWarpWeight(double q2, std::string channel) const;

#endif
