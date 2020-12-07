//==============================================================================
// This script just tests how TargetUtils works.
// TU is used for target counting for the xsec normalization step
//==============================================================================
#include <iostream>

#include "PlotUtils/TargetUtils.h"

void NTargetsChecks() {
  const double apothem = 865.;  // cm
  const double upstream = 5900.;
  const double downstream = 8430.;
  const int nplanes =
      PlotUtils::TargetUtils::Get().GetNPlanes(upstream, downstream);
  std::cout << "N Planes in FidVol(" << upstream << "," << downstream
            << ") : " << nplanes << "\n";

  // Data and MC will be the same with NX processing:
  // "In NX the mc was generated with the as-built geometry"

  // Data
  std::cout << "Data\n";
  bool is_mc = false;
  double targets = PlotUtils::TargetUtils::Get().GetTrackerNNucleons(
      upstream, downstream, is_mc, apothem);
  double targets_check = PlotUtils::TargetUtils::Get().GetTrackerNNucleons(
      nplanes, is_mc, apothem);
  std::cout << "N Nucleons in FidVol(" << upstream << "," << downstream
            << ") : " << targets << "\n";
  std::cout << "N Nucleons in nplanes(" << nplanes
            << ") (cross check): " << targets_check << "\n";

  // MC
  std::cout << "MC\n";
  is_mc = true;
  targets = PlotUtils::TargetUtils::Get().GetTrackerNNucleons(
      upstream, downstream, is_mc, apothem);
  targets_check = PlotUtils::TargetUtils::Get().GetTrackerNNucleons(
      nplanes, is_mc, apothem);
  std::cout << "N Nucleons in FidVol(" << upstream << "," << downstream
            << ") : " << targets << "\n";
  std::cout << "N Nucleons in nplanes(" << nplanes
            << ") (cross check): " << targets_check << "\n";
}

// const double apothem = 865.;
// const double upstream = 5900.;
// const double downstream = 8430.;
// double targets = TargetUtils::Get().GetTrackerNNucleons(upstream, downstream,
// /*isMC =*/ false, apothem) int nplanes =
// TargetUtils::Get().GetNPlanes(upstream, downstream); double targets_check =
// TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, apothem);

// cout << "N Planes in FidVol(" << upstream << "," << downstream << ") : " <<
// nplanes << "\n"; cout << "N Nucleons in FidVol(" << upstream << "," <<
// downstream << ") : " << targets << "\n"; cout << "N Nucleons in FidVol(" <<
// upstream << "," << downstream << ") (cross check): " << targets_check << "\n";

// const int nplanes = 2 * ( 80 - 27 + 1 ); //fiducial volume -> modules 27-80
// double targets = TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/
// false, apothem); std::cout << "NPlanes 27-80: " << nplanes   << "\n";

// std::cout << "NNucleons 27-80: " << targets   << "\n";

// cout << "FidVol(" << upstream << "," << downstream << ") NNucleons: "
//     << TargetUtils::Get().GetTrackerNNucleons(upstream, downstream, /*isMC
//     =*/ false, apothem) << "\n";

// const int nplanes_2 = 2 * ( 79 - 27 + 1 ); //fiducial volume -> modules 27-79
// std::cout << "NPlanes 27-79: " << nplanes_2 << "\n";
// double targets_2 = TargetUtils::Get().GetTrackerNNucleons( nplanes_2, /*isMC
// =*/ false, apothem); std::cout << "NNucleons 27-79: " << targets_2 << "\n";

// cout << TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ true,
// apothem) << endl; cout << TargetUtils::Get().GetTrackerNProtons(     nplanes,
// /*isMC =*/ true, apothem) << endl; cout << "MC\t"            <<
// TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ true,  apothem) <<
// endl;//MC cout << "Data\t"          << TargetUtils::Get().GetTrackerNNucleons(
// nplanes, /*isMC =*/ false, apothem) << endl;//Data cout << "MC Neutrons\t" <<
// TargetUtils::Get().GetTrackerNNeutrons( nplanes, /*isMC =*/ true,  apothem) <<
// endl;//MC cout << "Data Neutrons\t" << TargetUtils::Get().GetTrackerNNeutrons(
// nplanes, /*isMC =*/ false, apothem) << endl;//Data
////Print out some useful information...
////double flux = flux_normalization->GetBinContent(1);
// double nC12 = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, true,
// apothem); cout << " number of targets: "      << targets      << " neutrons"
// << endl; cout << " number of C12(MC): "      << nC12         << endl; cout <<
// " ratio of nC12/neutrons: " << nC12/targets << endl;
////cout << " flux integrated: "        << flux         << endl;

// cout << TargetUtils::Get().GetTrackerNNeutrons(5980.,8422.,true,apothem) <<
// "\t"
//     << TargetUtils::Get().GetNPlanes(5980,8422)                      <<"'\t
//     This code." << endl;
