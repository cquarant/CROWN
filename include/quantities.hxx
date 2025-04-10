#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "basefunctions.hxx"
#include "defaults.hxx"
#include "utility/Logger.hxx"
#include "vectoroperations.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <Math/Boost.h>
// #include "../include/MT2/lester_mt2.hxx"

namespace quantities {
//mingxuan MT2
// ROOT::RDF::RNode calc_MT2(ROOT::RDF::RNode df, const std::string &outputname,
//                         const std::string &mu1_p4, const std::string &mu2_p4, const std::string &met_p4);
// //mingxuan add
ROOT::RDF::RNode calc_Mct(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &mu1_p4, const std::string &mu2_p4);
ROOT::RDF::RNode calc_ratio(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &A, const std::string &B);
ROOT::RDF::RNode calculateNeutrinoPz(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &lep_p4, const std::string &met_p4);    
ROOT::RDF::RNode calculate_kT(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p1, const std::string &p2);
ROOT::RDF::RNode calculate_antikT(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p1, const std::string &p2);
ROOT::RDF::RNode mT_MHT(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met);
ROOT::RDF::RNode deltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4);
ROOT::RDF::RNode deltaEta(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4);
ROOT::RDF::RNode pt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &inputvector);
ROOT::RDF::RNode eta(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector);
ROOT::RDF::RNode phi(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector);
ROOT::RDF::RNode mass(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &inputvector);
ROOT::RDF::RNode dxy(ROOT::RDF::RNode df, const std::string &outputname,
                     const int &position, const std::string &pairname,
                     const std::string &dxycolumn);
ROOT::RDF::RNode dz(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &dzcolumn);
ROOT::RDF::RNode charge(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &position, const std::string &pairname,
                        const std::string &chargecolumn);
ROOT::RDF::RNode scalarPtSum(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &pt_1, const std::string &pt_2, const std::string &pt_3);
ROOT::RDF::RNode scalarPtSum_met(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &pt_1, const std::string &pt_2, const std::string &pt_3 ,const std::string &pt_4 );
ROOT::RDF::RNode m_vis(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pt_W(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pt_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode deltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4);
ROOT::RDF::RNode deltaPhi_WH(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4, const std::string &p_3_p4);
// ROOT::RDF::RNode
// p4_fastmtt(ROOT::RDF::RNode df, const std::string &outputname,
//            const std::string &pt_1, const std::string &pt_2,
//            const std::string &eta_1, const std::string &eta_2,
//            const std::string &phi_1, const std::string &phi_2,
//            const std::string &mass_1, const std::string &mass_2,
//            const std::string &met_pt, const std::string &met_phi,
//            const std::string &met_cov_xx, const std::string &met_cov_xy,
//            const std::string &met_cov_yy, const std::string &decay_mode_1,
//            const std::string &decay_mode_2, const std::string &finalstate);
ROOT::RDF::RNode pt_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pzetamissvis(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &p_1_p4,
                              const std::string &p_2_p4,
                              const std::string &met);
ROOT::RDF::RNode mTdileptonMET(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &p_1_p4,
                               const std::string &p_2_p4,
                               const std::string &met);
ROOT::RDF::RNode deltaR(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4);
ROOT::RDF::RNode mT(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met);
ROOT::RDF::RNode pt_tt(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &p_1_p4, const std::string &p_2_p4,
                       const std::string &met);
ROOT::RDF::RNode pt_ttjj(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::string &p_1_p4, const std::string &p_2_p4,
                         const std::string &jet_1_p4,
                         const std::string &jet_2_p4, const std::string &met);
ROOT::RDF::RNode pt_dijet(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &jet_1_p4,
                          const std::string &jet_2_p4);
ROOT::RDF::RNode jet_hemisphere(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &jet_1_p4,
                                const std::string &jet_2_p4);
ROOT::RDF::RNode mt_tot(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met);
ROOT::RDF::RNode isolation(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &isolationcolumn);
ROOT::RDF::RNode pdgid(ROOT::RDF::RNode df, const std::string &outputname,
                       const int &position, const std::string &pairname,
                       const std::string &pdgidcolumn);
ROOT::RDF::RNode NumberOfGoodObjects(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const std::string &goodmask);
ROOT::RDF::RNode NumberOfGoodLeptons(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const std::string &goodleptons);
namespace tau {
ROOT::RDF::RNode decaymode(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &decaymodecolumn);
ROOT::RDF::RNode genmatch(ROOT::RDF::RNode df, const std::string &outputname,
                          const int &position, const std::string &pairname,
                          const std::string &genmatchcolumn);
ROOT::RDF::RNode matching_jet_pt(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const int &position,
                                 const std::string &pairname,
                                 const std::string &taujet_index,
                                 const std::string &jetpt_column);
ROOT::RDF::RNode matching_genjet_pt(
    ROOT::RDF::RNode df, const std::string &outputname, const int &position,
    const std::string &pairname, const std::string &taujet_index,
    const std::string &genjet_index, const std::string &genjetpt_column);
ROOT::RDF::RNode TauIDFlag(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &nameID, const int &idxID);
} // end namespace tau
namespace muon {
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn);
ROOT::RDF::RNode is_global(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &globalflagcolumn);
} // namespace muon
namespace electron {
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn);
} // end namespace electron
namespace boostedbbtt {
    ROOT::RDF::RNode dR_tau_fatjet0(
        ROOT::RDF::RNode df, 
        const std::string &dR_tau_fatjet,
        const std::string &p0_fatjet,   
        const std::string &taupt,     
        const std::string &taueta,   
        const std::string &tauphi, 
        const std::string &taumass);

    void DecomposeMomentum(
        const ROOT::Math::PtEtaPhiMVector& p, 
        const ROOT::Math::PtEtaPhiMVector& a, 
        const ROOT::Math::PtEtaPhiMVector& b,
        ROOT::Math::PtEtaPhiMVector& pa,
        ROOT::Math::PtEtaPhiMVector& pb);

    ROOT::RDF::RNode Mass_CA(
        ROOT::RDF::RNode df,
        const std::string &p0_fatjet,
        const std::string &p1_fatjet,
        const std::string &met_pt,
        const std::string &met_phi,
        const std::string &mass_caca);

    ROOT::RDF::RNode CA_ttMAss(
        ROOT::RDF::RNode df,
        const std::string &x0x1,
        const std::string &FatJet_mass,
        const std::string &result);
    
    ROOT::RDF::RNode CA_ttMAss_fatjet(
        ROOT::RDF::RNode df,
        const std::string &p0_fatjet,
        const std::string &met_pt,
        const std::string &met_phi,
        const std::string &result); 
    
    ROOT::RDF::RNode Score_9X(
        ROOT::RDF::RNode df,
        const std::string &score,
        const std::string &result);

    ROOT::RDF::RNode Score_original(
        ROOT::RDF::RNode df,
        const std::string &Score,
        const std::string &QCD,
        const std::string &result); 
         
} // end namespace boostedbbtt
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */
