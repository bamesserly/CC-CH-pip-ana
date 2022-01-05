
std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    delta_t = [](const CVUniverse& univ, const MichelEvent& evt,
                 const int whichMichel) {
      int evttype = evt.eventtype;
      double micheltime = evt.m_nmichels[whichMichel].time;
      double vtxtime = univ.GetVertex().t();
      double deltat =
          (micheltime -
           vtxtime / 1000.);  // hopefully this is in microseconds (mus)
      // if (evttype == 1) return deltat;
      // else return -9999.;
      return deltat;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    permichel_range = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      // if (evt.eventtype == 1) return micheldist;
      // else return -9999.;
      return micheldist;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    pertruepimichel_range = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (evt.m_nmichels[whichMichel].true_parentpdg == 211)
        return micheldist;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    permichel_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                       const int whichMichel) {
      return evt.m_nmichels[whichMichel].pionKE;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    overlay_vtx_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      double KE = evt.m_nmichels[whichMichel].pionKE;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      if (overlayfrac > .5 && (matchtype == 1 || matchtype == 2))
        return KE;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    overlay_clus_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                          const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      double KE = evt.m_nmichels[whichMichel].pionKE;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      if (overlayfrac > .5 && (matchtype == 3 || matchtype == 4))
        return KE;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_goodvtx_tpi = [](const CVUniverse& univ,
                                const MichelEvent& evt,
                                const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      double KE = evt.m_nmichels[whichMichel].pionKE;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd == recoEnd &&
          (matchtype == 1 || matchtype == 2))
        return KE;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_goodclus_tpi = [](const CVUniverse& univ,
                                 const MichelEvent& evt,
                                 const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      double KE = evt.m_nmichels[whichMichel].pionKE;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd == recoEnd &&
          (matchtype == 3 || matchtype == 4))
        return KE;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_badvtx_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      double KE = evt.m_nmichels[whichMichel].pionKE;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd != recoEnd &&
          (matchtype == 1 || matchtype == 2)) {
        // univ.PrintArachneLink();
        // std::cout << "Printing Michel Time for bad VTX match type "  <<
        // evt.m_nmichels[whichMichel].time << std::endl;
        return KE;
      } else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_badclus_tpi = [](const CVUniverse& univ,
                                const MichelEvent& evt,
                                const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      double KE = evt.m_nmichels[whichMichel].pionKE;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd != recoEnd &&
          (matchtype == 3 || matchtype == 4)) {
        // univ.PrintArachneLink();
        // std::cout << "Printing Michel Time for bad CLUS match type "  <<
        // evt.m_nmichels[whichMichel].time << std::endl;
        return KE;
      } else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    overlay_vtx_range = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;

      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      if (overlayfrac > .5 && (matchtype == 1 || matchtype == 2))
        return micheldist;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    overlay_clus_range = [](const CVUniverse& univ, const MichelEvent& evt,
                            const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;

      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      if (overlayfrac > .5 && (matchtype == 3 || matchtype == 4))
        return micheldist;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_goodvtx_range = [](const CVUniverse& univ,
                                  const MichelEvent& evt,
                                  const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;

      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd == recoEnd &&
          (matchtype == 1 || matchtype == 2))
        return micheldist;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> truemichel_goodclus_range 
                              = [](const CVUniverse& univ,
                                   const MichelEvent& evt,
                                   const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd == recoEnd &&
          (matchtype == 3 || matchtype == 4))
        return micheldist;
      else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_badvtx_range = [](const CVUniverse& univ,
                                 const MichelEvent& evt,
                                 const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd != recoEnd &&
          (matchtype == 1 || matchtype == 2)) {
        //      univ.PrintArachneLink();
        //	std::cout << "Printing Michel Time for bad VTX match type "  <<
        //evt.m_nmichels[whichMichel].time << std::endl;
        return micheldist;
      } else
        return -9999.;
    };
std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    truemichel_badclus_range = [](const CVUniverse& univ,
                                  const MichelEvent& evt,
                                  const int whichMichel) {
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
      int matchtype = evt.m_nmichels[whichMichel].BestMatch;
      int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
      int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
      if (overlayfrac < .5 && trueEnd != recoEnd &&
          (matchtype == 3 || matchtype == 4)) {
        //	 univ.PrintArachneLink();
        //       std::cout << "Printing Michel Time for bad CLUS match type "
        //       << evt.m_nmichels[whichMichel].time << std::endl;
        return micheldist;
      } else
        return -9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_XZ = [](const CVUniverse& univ, const MichelEvent& evt,
                   const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].best_XZ;
      return twoDdist;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_UZ = [](const CVUniverse& univ, const MichelEvent& evt,
                   const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].best_UZ;
      return twoDdist;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_VZ = [](const CVUniverse& univ, const MichelEvent& evt,
                   const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].best_VZ;
      return twoDdist;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_XZ_upvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 1)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_UZ_upvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 1)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_VZ_upvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 1)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_XZ_upclus = [](const CVUniverse& univ, const MichelEvent& evt,
                          const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 1)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_UZ_upclus = [](const CVUniverse& univ, const MichelEvent& evt,
                          const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 1)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_VZ_upclus = [](const CVUniverse& univ, const MichelEvent& evt,
                          const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 1)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_XZ_downclus = [](const CVUniverse& univ, const MichelEvent& evt,
                            const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 2)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_UZ_downclus = [](const CVUniverse& univ, const MichelEvent& evt,
                            const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 2)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_VZ_downclus = [](const CVUniverse& univ, const MichelEvent& evt,
                            const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 2)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_XZ_downvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 2)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_UZ_downvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 2)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    michel_VZ_downvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double twoDdist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
      int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
      if (trueend == 2)
        return twoDdist;
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    pion_angle = [](const CVUniverse& univ, const MichelEvent& evt,
                    const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].best_angle;
      return cos(angle);
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    pion_angle_range1 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].best_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist <= 150.)
        return cos(angle);
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    pion_angle_range2 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].best_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist > 150. && micheldist <= 250.)
        return cos(angle);
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    pion_angle_range3 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].best_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist > 250. && micheldist <= 500.)
        return cos(angle);
      else
        return 9999.;
    };
std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    pion_angle_range4 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].best_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist > 500.)
        return cos(angle);
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    true_angle = [](const CVUniverse& univ, const MichelEvent& evt,
                    const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].true_angle;
      return cos(angle);
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    true_angle_range1 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].true_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist <= 150.)
        return cos(angle);
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    true_angle_range2 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].true_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist > 150 && micheldist <= 250.)
        return cos(angle);
      else
        return 9999.;
    };
std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    true_angle_range3 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].true_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist > 250 && micheldist <= 500.)
        return cos(angle);
      else
        return 9999.;
    };

std::function<double(const CVUniverse&, const MichelEvent&, const int)>
    true_angle_range4 = [](const CVUniverse& univ, const MichelEvent& evt,
                           const int whichMichel) {
      double angle = evt.m_nmichels[whichMichel].true_angle;
      double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
      if (micheldist > 500)
        return cos(angle);
      else
        return 9999.;
    };
