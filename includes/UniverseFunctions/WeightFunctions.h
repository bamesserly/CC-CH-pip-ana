#ifndef WeightFunctions_h
#define WeightFunctions_h

// Get Weight
virtual double GetWeight() const;
virtual double GetDiffractiveWeight() const;
// Warping
virtual double GetGenieWarpWeight() const;
virtual double GetLowQ2PiWarpWeight(double q2, std::string channel) const;
virtual double GetAnisoDeltaDecayWarpWeight() const;

#endif
