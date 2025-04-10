#ifndef GUARDJETS_H
#define GUARDJETS_H

#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <typeinfo>

namespace jet {

// vhmm extend to N particle overlap removal
/// Function to veto jets overlapping with particle candidates
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask \param[in] jet_eta name
/// of the jet etas \param[in] jet_phi name of the jet phis \param[in] p4_1 four
/// vector of the first particle candidate \param[in] p4_2 four vector of the
/// second particle candidate \param[in] deltaRmin minimum required distance in
/// dR between jets and particle candidates
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &output_col,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &muon_eta, const std::string &muon_phi, const std::string &muon_coll,
                    const float &deltaRmin) {
    auto df1 = df.Define(
        output_col,
        [deltaRmin](const ROOT::RVec<float> &jet_eta,
                    const ROOT::RVec<float> &jet_phi,
                    const ROOT::RVec<float> &muon_eta,
                    const ROOT::RVec<float> &muon_phi,
                    const ROOT::RVec<int> &muon_coll) {
            Logger::get("VetoOverlappingJets (N particles)")
                ->debug("Checking jets");
            ROOT::RVec<int> mask(jet_eta.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_eta.at(idx),
                                                 jet_phi.at(idx));
                Logger::get("VetoOverlappingJets (N particles)")
                    ->debug("Jet {}:  Eta: {} Phi: {} ", idx, jet.Eta(), jet.Phi());
                int _mdx = 0;
                for(std::size_t mdx = 0; mdx < muon_coll.size(); ++mdx){
                    // input collection, not mask, so they are index, no need to judge mask 0 or 1
                    // need to notice that the muon index is muon_coll[mdx]
                    // if( muon_mask[mdx] == 0 ) continue; // only check with the selected muons
                    ROOT::Math::RhoEtaPhiVectorF muon(0, muon_eta.at(muon_coll[mdx]), muon_phi.at(muon_coll[mdx]));
                    Logger::get("VetoOverlappingJets (N particles)")
                        ->debug("Lepton {}:  Eta: {} Phi: {} ", _mdx, muon.Eta(), muon.Phi());
                    auto deltaR = ROOT::Math::VectorUtil::DeltaR(jet, muon);
                    Logger::get("VetoOverlappingJets (N particles)")
                        ->debug("DeltaR {}", deltaR);
                    mask[idx] = mask[idx]&&(deltaR > deltaRmin);
                    ++_mdx;
                }
            }
            Logger::get("VetoOverlappingJets (N particles)")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, muon_eta, muon_phi, muon_coll});
    return df1;
}



/// Function to veto jets overlapping with particle candidates
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask \param[in] jet_eta name
/// of the jet etas \param[in] jet_phi name of the jet phis \param[in] p4_1 four
/// vector of the first particle candidate \param[in] p4_2 four vector of the
/// second particle candidate \param[in] deltaRmin minimum required distance in
/// dR between jets and particle candidates
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &output_col,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &p4_1, const std::string &p4_2,
                    const float &deltaRmin) {
    auto df1 = df.Define(
        output_col,
        [deltaRmin](const ROOT::RVec<float> &jet_eta,
                    const ROOT::RVec<float> &jet_phi,
                    const ROOT::Math::PtEtaPhiMVector &p4_1,
                    const ROOT::Math::PtEtaPhiMVector &p4_2) {
            Logger::get("VetoOverlappingJets (2 particles)")
                ->debug("Checking jets");
            ROOT::RVec<int> mask(jet_eta.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_eta.at(idx),
                                                 jet_phi.at(idx));
                Logger::get("VetoOverlappingJets (2 particles)")
                    ->debug("Jet:  Eta: {} Phi: {} ", jet.Eta(), jet.Phi());
                Logger::get("VetoOverlappingJets (2 particles)")
                    ->debug("Letpon 1 {}:  Eta: {} Phi: {}, Pt{}", p4_1,
                            p4_1.Eta(), p4_1.Phi(), p4_1.Pt());
                Logger::get("VetoOverlappingJets (2 particles)")
                    ->debug("Lepton 2 {}:  Eta: {} Phi: {}, Pt{}", p4_2,
                            p4_2.Eta(), p4_2.Phi(), p4_2.Pt());
                auto deltaR_1 = ROOT::Math::VectorUtil::DeltaR(jet, p4_1);
                auto deltaR_2 = ROOT::Math::VectorUtil::DeltaR(jet, p4_2);
                Logger::get("VetoOverlappingJets (2 particles)")
                    ->debug("DeltaR 1 {}", deltaR_1);
                Logger::get("VetoOverlappingJets (2 particles)")
                    ->debug("DeltaR 2 {}", deltaR_2);
                mask[idx] = (deltaR_1 > deltaRmin && deltaR_2 > deltaRmin);
            }
            Logger::get("VetoOverlappingJets (2 particles)")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, p4_1, p4_2});
    return df1;
}

/// Function to veto jets overlapping with particle candidates
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask \param[in] jet_eta name
/// of the jet etas \param[in] jet_phi name of the jet phis \param[in] p4_1 four
/// vector of the first particle candidate \param[in] deltaRmin minimum required
/// distance in dR between jets and particle candidates
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &output_col,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &p4_1, const float &deltaRmin) {
    auto df1 = df.Define(
        output_col,
        [deltaRmin](const ROOT::RVec<float> &jet_eta,
                    const ROOT::RVec<float> &jet_phi,
                    const ROOT::Math::PtEtaPhiMVector &p4_1) {
            Logger::get("VetoOverlappingJets")->debug("Checking jets");
            ROOT::RVec<int> mask(jet_eta.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_eta.at(idx),
                                                 jet_phi.at(idx));
                Logger::get("VetoOverlappingJets")
                    ->debug("Jet:  Eta: {} Phi: {} ", jet.Eta(), jet.Phi());
                Logger::get("VetoOverlappingJets")
                    ->debug("Letpon 1 {}:  Eta: {} Phi: {}, Pt{}", p4_1,
                            p4_1.Eta(), p4_1.Phi(), p4_1.Pt());
                auto deltaR_1 = ROOT::Math::VectorUtil::DeltaR(jet, p4_1);
                Logger::get("VetoOverlappingJets")
                    ->debug("DeltaR 1 {}", deltaR_1);
                mask[idx] = (deltaR_1 > deltaRmin);
            }
            Logger::get("VetoOverlappingJets")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, p4_1});
    return df1;
}

/// Function to veto jets overlapping with particle candidates (with isolation
/// condition)
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask \param[in] jet_eta name
/// of the jet etas \param[in] jet_phi name of the jet phis \param[in] p4_1 four
/// vector of the first particle candidate \param[in] deltaRmin minimum required
/// distance in dR between jets and particle candidates
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode VetoOverlappingJetsIsoLepOnly(ROOT::RDF::RNode df,
                                               const std::string &output_col,
                                               const std::string &jet_eta,
                                               const std::string &jet_phi,
                                               const std::string &p4_1,
                                               const std::string &lep_is_iso,
                                               const float &deltaRmin) {
    auto df1 = df.Define(
        output_col,
        [deltaRmin](
            const ROOT::RVec<float> &jet_eta, const ROOT::RVec<float> &jet_phi,
            const ROOT::Math::PtEtaPhiMVector &p4_1, const int &lep_is_iso) {
            Logger::get("VetoOverlappingJets")->debug("Checking jets");
            ROOT::RVec<int> mask(jet_eta.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_eta.at(idx),
                                                 jet_phi.at(idx));
                Logger::get("VetoOverlappingJets")
                    ->debug("Jet:  Eta: {} Phi: {} ", jet.Eta(), jet.Phi());
                Logger::get("VetoOverlappingJets")
                    ->debug("Letpon 1 {}:  Eta: {} Phi: {}, Pt{}", p4_1,
                            p4_1.Eta(), p4_1.Phi(), p4_1.Pt());
                auto deltaR_1 = ROOT::Math::VectorUtil::DeltaR(jet, p4_1);
                Logger::get("VetoOverlappingJets")
                    ->debug("DeltaR 1 {}", deltaR_1);
                if (lep_is_iso == +1)
                    mask[idx] = (deltaR_1 > deltaRmin);
            }
            Logger::get("VetoOverlappingJets")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, p4_1, lep_is_iso});
    return df1;
}

/// Function to determine pt order of jets
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask
/// \param[in] jet_pt name of the jet pts
/// \param[in] jetmask_name name of the mask marking all valid
/// jets to be considered
///
/// \return a dataframe containing a list of jet indices sorted by pt
ROOT::RDF::RNode OrderJetsByPt(ROOT::RDF::RNode df,
                               const std::string &output_col,
                               const std::string &jet_pt,
                               const std::string &jetmask_name) {
    auto df1 = df.Define(
        output_col,
        [output_col, jetmask_name](const ROOT::RVec<int> &jetmask,
                                   const ROOT::RVec<float> &jet_pt) {
            Logger::get("OrderJetsByPt")
                ->debug("Ordering good jets from {} by pt, output stored in {}",
                        jetmask_name, output_col);
            Logger::get("OrderJetsByPt")->debug("Jetpt before {}", jet_pt);
            Logger::get("OrderJetsByPt")->debug("Mask {}", jetmask);
            auto good_jets_pt =
                ROOT::VecOps::Where(jetmask > 0, jet_pt, (float)0.);
            Logger::get("OrderJetsByPt")->debug("Jetpt after {}", good_jets_pt);
            // we have to convert the result into an RVec of ints since argsort
            // gives back an unsigned long vector
            auto temp = ROOT::VecOps::Intersect(
                ROOT::VecOps::Argsort(good_jets_pt,
                                      [](double x, double y) { return x > y; }),
                ROOT::VecOps::Nonzero(good_jets_pt));
            Logger::get("OrderJetsByPt")->debug("jet Indices {}", temp);
            ROOT::RVec<int> result(temp.size());
            std::transform(temp.begin(), temp.end(), result.begin(),
                           [](unsigned long int x) { return (int)x; });
            Logger::get("OrderJetsByPt")->debug("jet Indices int {}", result);
            return result;
        },
        {jetmask_name, jet_pt});
    return df1;
}

//jin for double sort
ROOT::RDF::RNode OrderJetsByCustomCriteria(ROOT::RDF::RNode df,
                                           const std::string &output_col,
                                           const std::string &jet_xtt,
                                           const std::string &jet_xbb,
                                           const std::string &jetmask_name) {
    auto df1 = df.Define(
        output_col,
        [output_col, jetmask_name](const ROOT::RVec<int> &jetmask,
                                   const ROOT::RVec<float> &jet_xtt,
                                   const ROOT::RVec<float> &jet_xbb) {
            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Ordering good jets by jet_xtt and jet_xbb, output stored in {}",
                        output_col);

            // 1. First, sort jets by jet_xtt (max jet_xtt first)
            auto good_jets_xtt = ROOT::VecOps::Where(jetmask > 0, jet_xtt, (float)0.);
            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Jetxtt before sorting: {}", good_jets_xtt);

            auto sorted_by_xtt_indices = ROOT::VecOps::Intersect(
                ROOT::VecOps::Argsort(good_jets_xtt, [](double x, double y) { return x > y; }),
                ROOT::VecOps::Nonzero(good_jets_xtt));

            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Sorted by jet_xtt indices: {}", sorted_by_xtt_indices);

            // 2. Get the index of the jet with the maximum jet_xtt
            int max_xtt_index = sorted_by_xtt_indices[0];
            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Max jet_xtt index: {}", max_xtt_index);

            // 3. Now remove the jet with max jet_xtt and sort the rest by jet_xbb (max jet_xbb second)
            ROOT::RVec<int> remaining_indices;
            for (size_t i = 1; i < sorted_by_xtt_indices.size(); ++i) {
                remaining_indices.push_back(sorted_by_xtt_indices[i]);
            }

            auto remaining_jets_xbb = ROOT::VecOps::Where(jetmask > 0, jet_xbb, (float)0.);
            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Jetxbb before sorting: {}", remaining_jets_xbb);

            auto sorted_by_xbb_indices = ROOT::VecOps::Intersect(
                ROOT::VecOps::Argsort(remaining_jets_xbb, [](double x, double y) { return x > y; }),
                ROOT::VecOps::Nonzero(remaining_jets_xbb));

            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Sorted by jet_xbb indices: {}", sorted_by_xbb_indices);

            // 4. Combine the indices: first the max jet_xtt, then sorted remaining by jet_xbb
            ROOT::RVec<int> final_order;
            final_order.push_back(max_xtt_index);
            final_order.insert(final_order.end(), sorted_by_xbb_indices.begin(), sorted_by_xbb_indices.end());

            Logger::get("OrderJetsByCustomCriteria")
                ->debug("Final jet order: {}", final_order);

            // Return the final ordered jet indices
            return final_order;
        },
        {jetmask_name, jet_xtt, jet_xbb});
    return df1;
}


} // end namespace jet

namespace physicsobject {
namespace jet {

/// Function to cut jets based on the jet ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID, const int &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterJetID(idxID), {nameID});
    return df1;
}
ROOT::RDF::RNode CutUCharID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID, const unsigned char &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterJetUCharID(idxID), {nameID});
    return df1;
}
/// Function to cut jets based on the jet pileup ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed \param[in] jet_pt name of the
/// input jet pts \param[in] jet_pt_cut threshold for the input jet pts
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutPUID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const std::string &jet_pt,
                         const int &idxID, const float &jet_pt_cut) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterJetPUID(idxID, jet_pt_cut),
                  {nameID, jet_pt});
    return df1;
}

/// function to do jet veto map for data (and mc)
ROOT::RDF::RNode
JetVetoMap(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                const std::string &jet_pt, const std::string &jet_eta,
                const std::string &jet_phi,
                const std::string &jet_veto_map, const std::string &jet_veto_tag) {
    // loading jet veto maps
    auto jet_veto_map_evaluator = correction::CorrectionSet::from_file(jet_veto_map)->at(
        jet_veto_tag);
    auto jet_veto_SF =
        [jet_veto_map_evaluator](const float eta, const float phi) {
            auto tmp_phi  = phi;
            if (phi > 3.141592653589793) tmp_phi = phi - (3.141592653589793 * 2);
            if (phi < -3.141592653589793) tmp_phi = phi + (3.141592653589793 * 2);
            if (std::abs(eta) <= 5.191) { // 5.191 is the bin dege of the veto map
                // Non-zero value for (eta, phi) indicates that the region is vetoed.
                // link: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
                return jet_veto_map_evaluator->evaluate({ "jetvetomap", eta,  tmp_phi});
            } else {
                // if the jet eta out of the bin edge, no need to do the veto
                return 0.0;
            }
    };

    // lambda run with dataframe
    auto JetEnergyCorrectionLambda = [jet_veto_SF](
                                         const ROOT::RVec<float> &pt_values,
                                         const ROOT::RVec<float> &eta_values,
                                         const ROOT::RVec<float> &phi_values) {
        ROOT::RVec<float> pt_values_corrected;
        // return original pt firstly
        pt_values_corrected = pt_values;
        // // apply jet veto map. If any jet lies within jet veto map, reject the events (jet?). 
        // // at the object level it's not straightforward to veto the events, so we return a RVec of pt -999 for all jets 
        float pt_veto = -999.0;
        // Loop to check if any non-zero jet_veto_sf_value exists
        for (int i = 0; i < pt_values.size(); i++) {

            Logger::get("JetEnergyResolution")
                ->debug("checking jet veto map for index {} ", i);
            float jet_veto_sf_value = jet_veto_SF(eta_values.at(i), phi_values.at(i));
            if (jet_veto_sf_value != 0) {
                // find a non-zero value, since the jet in the veto regions, enforce the jet pt == -999.0
                pt_values_corrected.at(i) = pt_veto;
            }
        }
        return pt_values_corrected;
    };
    auto df1 = df.Define(corrected_jet_pt, JetEnergyCorrectionLambda,
                         {jet_pt, jet_eta, jet_phi});
    return df1;
}

/// Function to shift and smear jet pt for MC
///
/// \param[in] df the input dataframe
/// \param[out] corrected_jet_pt the name of the shifted and smeared jet pts
/// \param[in] jet_pt name of the input jet pts
/// \param[in] jet_eta name of the jet etas
/// \param[in] jet_phi name of the jet phis
/// \param[in] jet_area name of the jet catchment area
/// \param[in] jet_rawFactor name of the raw factor for jet pt
/// \param[in] jet_ID name of the jet ID
/// \param[in] gen_jet_pt name of the gen jet pts
/// \param[in] gen_jet_eta name of the gen jet etas
/// \param[in] gen_jet_phi name of the gen jet phis
/// \param[in] rho name of the pileup density
/// \param[in] reapplyJES boolean for reapplying the JES correction
/// \param[in] jes_shift_sources vector of JEC unc source names to be applied
/// in one group
/// \param[in] jes_shift parameter to control jet energy
/// scale shift: 0 - nominal; 1 - Up; -1 - Down
/// \param[in] jer_shift parameter to control jet energy resolution
/// shift: "nom"; "up"; "down"
/// \param[in] jec_file path to the file with JES/JER information
/// \param[in] jer_tag era dependent tag for JER
/// \param[in] jes_tag era dependent tag for JES
/// \param[in] jec_algo algorithm used for jets e.g. AK4PFchs
///
/// \return a dataframe containing the modified jet pts
ROOT::RDF::RNode
JetPtCorrection(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                const std::string &jet_pt, const std::string &jet_eta,
                const std::string &jet_phi, const std::string &jet_area,
                const std::string &jet_rawFactor, const std::string &jet_ID,
                const std::string &gen_jet_pt, const std::string &gen_jet_eta,
                const std::string &gen_jet_phi, const std::string &rho,
                bool reapplyJES,
                const std::vector<std::string> &jes_shift_sources,
                const int &jes_shift, const std::string &jer_shift,
                const std::string &jec_file, const std::string &jer_tag,
                const std::string &jes_tag, const std::string &jec_algo,
                const std::string &jet_veto_map, const std::string &jet_veto_tag) {
    // identifying jet radius from algorithm
    float jet_dR = 0.4;
    if (jec_algo.find("AK8") != std::string::npos) {
        jet_dR = 0.8;
    }
    // loading JES variations
    std::vector<std::shared_ptr<const correction::Correction>>
        JetEnergyScaleShifts;
    for (const auto &source : jes_shift_sources) {
        // check if any JES shift is chosen
        if (source != "" && source != "HEMIssue") {
            auto JES_source_evaluator =
                correction::CorrectionSet::from_file(jec_file)->at(
                    jes_tag + "_" + source + "_" + jec_algo);
            JetEnergyScaleShifts.push_back(JES_source_evaluator);
        }
    };
    // loading jet energy correction scale factor evaluation function
    auto JES_evaluator =
        correction::CorrectionSet::from_file(jec_file)->compound().at(
            jes_tag + "_L1L2L3Res_" + jec_algo);
    auto JetEnergyScaleSF = [JES_evaluator](const float area, const float eta,
                                            const float pt, const float rho) {
        if (std::abs(eta) < 4.7) {
            return JES_evaluator->evaluate({area, eta, pt, rho});
        } else {
            return 1.0;
        }
    };
    // loading relative pT resolution evaluation function
    auto JER_resolution_evaluator =
        correction::CorrectionSet::from_file(jec_file)->at(
            jer_tag + "_PtResolution_" + jec_algo);
    auto JetEnergyResolution = [JER_resolution_evaluator](const float eta,
                                                          const float pt,
                                                          const float rho) {
        if (std::abs(eta) < 4.7) {
            return JER_resolution_evaluator->evaluate({eta, pt, rho});
        } else {
            return 1.0;
        }
    };
    // loading JER scale factor evaluation function
    auto JER_SF_evaluator = correction::CorrectionSet::from_file(jec_file)->at(
        jer_tag + "_ScaleFactor_" + jec_algo);
    auto JetEnergyResolutionSF =
        [JER_SF_evaluator](const float eta, float pt, const std::string jer_shift) {
            try{ 
                // annoyingly the 2022 EE has 3 inputs?
                if (std::abs(eta) < 4.7) {
                    return JER_SF_evaluator->evaluate({eta, jer_shift});
                } else {
                    return 1.0;
                }
            }
            catch (const std::exception&) {
            if (std::abs(eta) < 4.7) {
                return JER_SF_evaluator->evaluate({eta, pt, jer_shift});
            } else {
                return 1.0;
            }
        }
    };
    // loading jet veto maps
    auto jet_veto_map_evaluator = correction::CorrectionSet::from_file(jet_veto_map)->at(
        jet_veto_tag);
    auto jet_veto_SF =
        [jet_veto_map_evaluator](const float eta, const float phi) {
            auto tmp_phi  = phi;
            if (phi > 3.141592653589793) tmp_phi = phi - (3.141592653589793 * 2);
            if (phi < -3.141592653589793) tmp_phi = phi + (3.141592653589793 * 2);
            if (std::abs(eta) <= 5.191) { // 5.191 is the bin dege of the veto map
                // Non-zero value for (eta, phi) indicates that the region is vetoed.
                // link: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
                return jet_veto_map_evaluator->evaluate({ "jetvetomap", eta,  tmp_phi});
            } else {
                // if the jet eta out of the bin edge, no need to do the veto
                return 0.0;
            }
    };

    // lambda run with dataframe
    auto JetEnergyCorrectionLambda = [reapplyJES, JetEnergyScaleShifts,
                                      JetEnergyScaleSF, JetEnergyResolution,
                                      JetEnergyResolutionSF, jes_shift_sources,
                                      jes_shift, jer_shift, jet_dR, jet_veto_SF](
                                         const ROOT::RVec<float> &pt_values,
                                         const ROOT::RVec<float> &eta_values,
                                         const ROOT::RVec<float> &phi_values,
                                         const ROOT::RVec<float> &area_values,
                                         const ROOT::RVec<float>
                                             &rawFactor_values,
                                         const ROOT::RVec<unsigned char> &ID_values,
                                         const ROOT::RVec<float> &gen_pt_values,
                                         const ROOT::RVec<float>
                                             &gen_eta_values,
                                         const ROOT::RVec<float>
                                             &gen_phi_values,
                                         const float &rho_value) {
        // random value generator for jet smearing
        TRandom3 randm = TRandom3(12345);
        ROOT::RVec<float> pt_values_corrected;
        for (int i = 0; i < pt_values.size(); i++) {
            float corr_pt = pt_values.at(i);
            if (reapplyJES) {
                // reapplying the JES correction
                float raw_pt = pt_values.at(i) * (1 - rawFactor_values.at(i));
                float corr = JetEnergyScaleSF(
                    area_values.at(i), eta_values.at(i), raw_pt, rho_value);
                corr_pt = raw_pt * corr;
                Logger::get("JetEnergyScale")
                    ->debug("reapplying JE scale: orig. jet pt {} to raw "
                            "jet pt {} to recorr. jet pt {}",
                            pt_values.at(i), raw_pt, corr_pt);
            }
            pt_values_corrected.push_back(corr_pt);

            // apply jet energy smearing - hybrid method as described in
            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
            float reso = JetEnergyResolution(
                eta_values.at(i), pt_values_corrected.at(i), rho_value);
            float resoSF = JetEnergyResolutionSF(eta_values.at(i), pt_values.at(i), jer_shift);
            Logger::get("JetEnergyResolution")
                ->debug("Calculate JER {}:  SF: {} resolution: {} ", jer_shift,
                        resoSF, reso);
            // gen jet matching algorithm for JER
            ROOT::Math::RhoEtaPhiVectorF jet(
                pt_values_corrected.at(i), eta_values.at(i), phi_values.at(i));
            float genjetpt = -1.0;
            Logger::get("JetEnergyResolution")
                ->debug("Going to smear jet:  Eta: {} Phi: {} ", jet.Eta(),
                        jet.Phi());
            double min_dR = std::numeric_limits<double>::infinity();
            for (int j = 0; j < gen_pt_values.size(); j++) {
                ROOT::Math::RhoEtaPhiVectorF genjet(gen_pt_values.at(j),
                                                    gen_eta_values.at(j),
                                                    gen_phi_values.at(j));
                Logger::get("JetEnergyResolution")
                    ->debug("Checking gen Jet:  Eta: {} Phi: {}", genjet.Eta(),
                            genjet.Phi());
                auto deltaR = ROOT::Math::VectorUtil::DeltaR(jet, genjet);
                if (deltaR > min_dR)
                    continue;
                if (deltaR < (jet_dR / 2.) &&
                    std::abs(pt_values_corrected.at(i) - gen_pt_values.at(j)) <
                        (3.0 * reso * pt_values_corrected.at(i))) {
                    min_dR = deltaR;
                    genjetpt = gen_pt_values.at(j);
                }
            }
            // if jet matches a gen jet scaling method is applied,
            // otherwise stochastic method
            if (genjetpt > 0.0) {
                Logger::get("JetEnergyResolution")
                    ->debug("Found gen jet for hybrid smearing method");
                double shift = (resoSF - 1.0) *
                               (pt_values_corrected.at(i) - genjetpt) /
                               pt_values_corrected.at(i);
                pt_values_corrected.at(i) *= std::max(0.0, 1.0 + shift);
            } else {
                Logger::get("JetEnergyResolution")
                    ->debug("No gen jet found. Applying stochastic smearing.");
                double shift = randm.Gaus(0, reso) *
                               std::sqrt(std::max(resoSF * resoSF - 1., 0.0));
                // if need, turn off the random for shift
                pt_values_corrected.at(i) *= std::max(0.0, 1.0 + shift);
            }
            Logger::get("JetEnergyResolution")
                ->debug("Shifting jet pt from {} to {} ", corr_pt,
                        pt_values_corrected.at(i));

            // apply uncertainty shifts related to the jet energy scale
            // mostly following
            // https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
            float pt_scale_sf = 1.0;
            if (jes_shift != 0.0) {
                if (jes_shift_sources.at(0) != "HEMIssue") {
                    // Differentiate between single source and combined source
                    // for reduced scheme
                    if (JetEnergyScaleShifts.size() == 1) {
                        pt_scale_sf =
                            1. +
                            jes_shift * JetEnergyScaleShifts.at(0)->evaluate(
                                            {eta_values.at(i),
                                             pt_values_corrected.at(i)});
                        Logger::get("JetEnergyScaleShift")
                            ->debug("Shifting jet pt by {} for single source "
                                    "with SF {}",
                                    jes_shift, pt_scale_sf);
                    } else {
                        float quad_sum = 0.;
                        for (const auto &evaluator : JetEnergyScaleShifts) {
                            quad_sum +=
                                std::pow(evaluator->evaluate(
                                             {eta_values.at(i),
                                              pt_values_corrected.at(i)}),
                                         2.0);
                        }
                        pt_scale_sf = 1. + jes_shift * std::sqrt(quad_sum);
                        Logger::get("JetEnergyScaleShift")
                            ->debug("Shifting jet pt by {} for multiple "
                                    "sources with SF {}",
                                    jes_shift, pt_scale_sf);
                    }
                }
                // for reference:
                // https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
                else if (jes_shift_sources.at(0) == "HEMIssue") {
                    if (jes_shift == (-1.) && pt_values_corrected.at(i) > 15. &&
                        phi_values.at(i) > (-1.57) &&
                        phi_values.at(i) < (-0.87) && ID_values.at(i) == 2) {
                        if (eta_values.at(i) > (-2.5) &&
                            eta_values.at(i) < (-1.3))
                            pt_scale_sf = 0.8;
                        else if (eta_values.at(i) > (-3.) &&
                                 eta_values.at(i) <= (-2.5))
                            pt_scale_sf = 0.65;
                    }
                }
            }
            pt_values_corrected.at(i) *= pt_scale_sf;
            Logger::get("JetEnergyScaleShift")
                ->debug("Shifting jet pt from {} to {} ",
                        pt_values_corrected.at(i) / pt_scale_sf,
                        pt_values_corrected.at(i));

            // if (pt_values_corrected.at(i)>15.0), this
            // correction should be propagated to MET
            // (requirement for type I corrections)
        }
        // // apply jet veto map. If any jet lies within jet veto map, reject the events (jet?). 
        // // at the object level it's not straightforward to veto the events, so we return a RVec of pt -999 for all jets 
        float pt_veto = -999.0;
        // Loop to check if any non-zero jet_veto_sf_value exists
        for (int i = 0; i < pt_values.size(); i++) {

            Logger::get("JetEnergyResolution")
                ->debug("checking jet veto map for index {} ", i);
            float jet_veto_sf_value = jet_veto_SF(eta_values.at(i), phi_values.at(i));
            if (jet_veto_sf_value != 0) {
                // find a non-zero value, since the jet in the veto regions, enforce the jet pt == -999.0
                pt_values_corrected.at(i) = pt_veto;
            }
        }
        return pt_values_corrected;
    };
    auto df1 = df.Define(corrected_jet_pt, JetEnergyCorrectionLambda,
                         {jet_pt, jet_eta, jet_phi, jet_area, jet_rawFactor,
                          jet_ID, gen_jet_pt, gen_jet_eta, gen_jet_phi, rho});
    return df1;
}

ROOT::RDF::RNode
JetPtCorrection_run2(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                const std::string &jet_pt, const std::string &jet_eta,
                const std::string &jet_phi, const std::string &jet_area,
                const std::string &jet_rawFactor, const std::string &jet_ID,
                const std::string &gen_jet_pt, const std::string &gen_jet_eta,
                const std::string &gen_jet_phi, const std::string &rho,
                bool reapplyJES,
                const std::vector<std::string> &jes_shift_sources,
                const int &jes_shift, const std::string &jer_shift,
                const std::string &jec_file, const std::string &jer_tag,
                const std::string &jes_tag, const std::string &jec_algo,
                const std::string &jet_veto_map, const std::string &jet_veto_tag) {
    // identifying jet radius from algorithm
    float jet_dR = 0.4;
    if (jec_algo.find("AK8") != std::string::npos) {
        jet_dR = 0.8;
    }
    // loading JES variations
    std::vector<std::shared_ptr<const correction::Correction>>
        JetEnergyScaleShifts;
    for (const auto &source : jes_shift_sources) {
        // check if any JES shift is chosen
        if (source != "" && source != "HEMIssue") {
            auto JES_source_evaluator =
                correction::CorrectionSet::from_file(jec_file)->at(
                    jes_tag + "_" + source + "_" + jec_algo);
            JetEnergyScaleShifts.push_back(JES_source_evaluator);
        }
    };
    // loading jet energy correction scale factor evaluation function
    auto JES_evaluator =
        correction::CorrectionSet::from_file(jec_file)->compound().at(
            jes_tag + "_L1L2L3Res_" + jec_algo);
    auto JetEnergyScaleSF = [JES_evaluator](const float area, const float eta,
                                            const float pt, const float rho) {
        if (std::abs(eta) < 4.7) {
            return JES_evaluator->evaluate({area, eta, pt, rho});
        } else {
            return 1.0;
        }
    };
    // loading relative pT resolution evaluation function
    auto JER_resolution_evaluator =
        correction::CorrectionSet::from_file(jec_file)->at(
            jer_tag + "_PtResolution_" + jec_algo);
    auto JetEnergyResolution = [JER_resolution_evaluator](const float eta,
                                                          const float pt,
                                                          const float rho) {
        if (std::abs(eta) < 4.7) {
            return JER_resolution_evaluator->evaluate({eta, pt, rho});
        } else {
            return 1.0;
        }
    };
    // loading JER scale factor evaluation function
    auto JER_SF_evaluator = correction::CorrectionSet::from_file(jec_file)->at(
        jer_tag + "_ScaleFactor_" + jec_algo);
    auto JetEnergyResolutionSF =
        [JER_SF_evaluator](const float eta, const std::string jer_shift) {
            if (std::abs(eta) < 4.7) {
                return JER_SF_evaluator->evaluate({eta, jer_shift});
            } else {
                return 1.0;
            }
    };
    // loading jet veto maps
    auto jet_veto_map_evaluator = correction::CorrectionSet::from_file(jet_veto_map)->at(
        jet_veto_tag);
    auto jet_veto_SF =
        [jet_veto_map_evaluator](const float eta, const float phi) {
            auto tmp_phi  = phi;
            if (phi > 3.141592653589793) tmp_phi = phi - (3.141592653589793 * 2);
            if (phi < -3.141592653589793) tmp_phi = phi + (3.141592653589793 * 2);
            if (std::abs(eta) <= 5.191) { // 5.191 is the bin dege of the veto map
                // Non-zero value for (eta, phi) indicates that the region is vetoed.
                // link: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
                return jet_veto_map_evaluator->evaluate({ "jetvetomap", eta,  tmp_phi});
            } else {
                // if the jet eta out of the bin edge, no need to do the veto
                return 0.0;
            }
    };

    // lambda run with dataframe
    auto JetEnergyCorrectionLambda_run2 = [reapplyJES, JetEnergyScaleShifts,
                                      JetEnergyScaleSF, JetEnergyResolution,
                                      JetEnergyResolutionSF, jes_shift_sources,
                                      jes_shift, jer_shift, jet_dR, jet_veto_SF](
                                         const ROOT::RVec<float> &pt_values,
                                         const ROOT::RVec<float> &eta_values,
                                         const ROOT::RVec<float> &phi_values,
                                         const ROOT::RVec<float> &area_values,
                                         const ROOT::RVec<float>
                                             &rawFactor_values,
                                         const ROOT::RVec<int> &ID_values,
                                         const ROOT::RVec<float> &gen_pt_values,
                                         const ROOT::RVec<float>
                                             &gen_eta_values,
                                         const ROOT::RVec<float>
                                             &gen_phi_values,
                                         const float &rho_value) {
        // random value generator for jet smearing
        TRandom3 randm = TRandom3(12345);
        ROOT::RVec<float> pt_values_corrected;        
        for (int i = 0; i < pt_values.size(); i++) {
            float corr_pt = pt_values.at(i);
            if (reapplyJES) {
                // reapplying the JES correction
                float raw_pt = pt_values.at(i) * (1 - rawFactor_values.at(i));
                float corr = JetEnergyScaleSF(
                    area_values.at(i), eta_values.at(i), raw_pt, rho_value);
                corr_pt = raw_pt * corr;
                Logger::get("JetEnergyScale")
                    ->debug("reapplying JE scale: orig. jet pt {} to raw "
                            "jet pt {} to recorr. jet pt {}",
                            pt_values.at(i), raw_pt, corr_pt);
            }
            pt_values_corrected.push_back(corr_pt);

            // apply jet energy smearing - hybrid method as described in
            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
            float reso = JetEnergyResolution(
                eta_values.at(i), pt_values_corrected.at(i), rho_value);
            // float resoSF = JetEnergyResolutionSF(eta_values.at(i), pt_values.at(i), jer_shift);
            float resoSF = JetEnergyResolutionSF(eta_values.at(i), jer_shift);
            Logger::get("JetEnergyResolution")
                ->debug("Calculate JER {}:  SF: {} resolution: {} ", jer_shift,
                        resoSF, reso);
            // gen jet matching algorithm for JER
            ROOT::Math::RhoEtaPhiVectorF jet(
                pt_values_corrected.at(i), eta_values.at(i), phi_values.at(i));
            float genjetpt = -1.0;
            Logger::get("JetEnergyResolution")
                ->debug("Going to smear jet:  Eta: {} Phi: {} ", jet.Eta(),
                        jet.Phi());
            double min_dR = std::numeric_limits<double>::infinity();
            for (int j = 0; j < gen_pt_values.size(); j++) {
                ROOT::Math::RhoEtaPhiVectorF genjet(gen_pt_values.at(j),
                                                    gen_eta_values.at(j),
                                                    gen_phi_values.at(j));
                Logger::get("JetEnergyResolution")
                    ->debug("Checking gen Jet:  Eta: {} Phi: {}", genjet.Eta(),
                            genjet.Phi());
                auto deltaR = ROOT::Math::VectorUtil::DeltaR(jet, genjet);
                if (deltaR > min_dR)
                    continue;
                if (deltaR < (jet_dR / 2.) &&
                    std::abs(pt_values_corrected.at(i) - gen_pt_values.at(j)) <
                        (3.0 * reso * pt_values_corrected.at(i))) {
                    min_dR = deltaR;
                    genjetpt = gen_pt_values.at(j);
                }
            }
            // if jet matches a gen jet scaling method is applied,
            // otherwise stochastic method
            if (genjetpt > 0.0) {
                Logger::get("JetEnergyResolution")
                    ->debug("Found gen jet for hybrid smearing method");
                double shift = (resoSF - 1.0) *
                               (pt_values_corrected.at(i) - genjetpt) /
                               pt_values_corrected.at(i);
                pt_values_corrected.at(i) *= std::max(0.0, 1.0 + shift);
            } else {
                Logger::get("JetEnergyResolution")
                    ->debug("No gen jet found. Applying stochastic smearing.");
                double shift = randm.Gaus(0, reso) *
                               std::sqrt(std::max(resoSF * resoSF - 1., 0.0));
                // if need, turn off the random for shift
                pt_values_corrected.at(i) *= std::max(0.0, 1.0 + shift);
            }
            Logger::get("JetEnergyResolution")
                ->debug("Shifting jet pt from {} to {} ", corr_pt,
                        pt_values_corrected.at(i));

            // apply uncertainty shifts related to the jet energy scale
            // mostly following
            // https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
            float pt_scale_sf = 1.0;
            if (jes_shift != 0.0) {
                if (jes_shift_sources.at(0) != "HEMIssue") {
                    // Differentiate between single source and combined source
                    // for reduced scheme
                    if (JetEnergyScaleShifts.size() == 1) {
                        pt_scale_sf =
                            1. +
                            jes_shift * JetEnergyScaleShifts.at(0)->evaluate(
                                            {eta_values.at(i),
                                             pt_values_corrected.at(i)});
                        Logger::get("JetEnergyScaleShift")
                            ->debug("Shifting jet pt by {} for single source "
                                    "with SF {}",
                                    jes_shift, pt_scale_sf);
                    } else {
                        float quad_sum = 0.;
                        for (const auto &evaluator : JetEnergyScaleShifts) {
                            quad_sum +=
                                std::pow(evaluator->evaluate(
                                             {eta_values.at(i),
                                              pt_values_corrected.at(i)}),
                                         2.0);
                        }
                        pt_scale_sf = 1. + jes_shift * std::sqrt(quad_sum);
                        Logger::get("JetEnergyScaleShift")
                            ->debug("Shifting jet pt by {} for multiple "
                                    "sources with SF {}",
                                    jes_shift, pt_scale_sf);
                    }
                }
                // for reference:
                // https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
                else if (jes_shift_sources.at(0) == "HEMIssue") {
                    if (jes_shift == (-1.) && pt_values_corrected.at(i) > 15. &&
                        phi_values.at(i) > (-1.57) &&
                        phi_values.at(i) < (-0.87) && ID_values.at(i) == 2) {
                        if (eta_values.at(i) > (-2.5) &&
                            eta_values.at(i) < (-1.3))
                            pt_scale_sf = 0.8;
                        else if (eta_values.at(i) > (-3.) &&
                                 eta_values.at(i) <= (-2.5))
                            pt_scale_sf = 0.65;
                    }
                }
            }
            pt_values_corrected.at(i) *= pt_scale_sf;
            Logger::get("JetEnergyScaleShift")
                ->debug("Shifting jet pt from {} to {} ",
                        pt_values_corrected.at(i) / pt_scale_sf,
                        pt_values_corrected.at(i));

            // if (pt_values_corrected.at(i)>15.0), this
            // correction should be propagated to MET
            // (requirement for type I corrections)
        }
        // // apply jet veto map. If any jet lies within jet veto map, reject the events (jet?). 
        // // at the object level it's not straightforward to veto the events, so we return a RVec of pt -999 for all jets 
        float pt_veto = -999.0;
        // Loop to check if any non-zero jet_veto_sf_value exists
        for (int i = 0; i < pt_values.size(); i++) {

            Logger::get("JetEnergyResolution")
                ->debug("checking jet veto map for index {} ", i);
            float jet_veto_sf_value = jet_veto_SF(eta_values.at(i), phi_values.at(i));
            if (jet_veto_sf_value != 0) {
                // find a non-zero value, since the jet in the veto regions, enforce the jet pt == -999.0
                pt_values_corrected.at(i) = pt_veto;
            }
        }
        return pt_values_corrected;
    };
    auto df1 = df.Define(corrected_jet_pt, JetEnergyCorrectionLambda_run2,
                         {jet_pt, jet_eta, jet_phi, jet_area, jet_rawFactor,
                          jet_ID, gen_jet_pt, gen_jet_eta, gen_jet_phi, rho});
    return df1;
}
/// Function to correct jet energy for data
///
/// \param[in] df the input dataframe
/// \param[out] corrected_jet_pt the name of the corrected jet pts
/// \param[in] jet_pt name of the input jet pts
/// \param[in] jet_eta name of the jet etas
/// \param[in] jet_area name of the jet catchment area
/// \param[in] jet_rawFactor name of the raw factor for jet pt
/// \param[in] rho name of the pileup density
/// \param[in] jec_file path to the file with jet energy correction information
/// \param[in] jes_tag era dependent tag for JES correction
/// \param[in] jec_algo algorithm used for jets e.g. AK4PFchs
///
/// \return a dataframe containing the modified jet pts
ROOT::RDF::RNode
JetPtCorrection_data(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                     const std::string &jet_pt, const std::string &jet_eta,
                     const std::string &jet_area,
                     const std::string &jet_rawFactor, const std::string &rho,
                     const std::string &jec_file, const std::string &jes_tag,
                     const std::string &jec_algo) {
    if (jes_tag != "") {
        // loading jet energy correction scale factor evaluation function
        auto JES_evaluator =
            correction::CorrectionSet::from_file(jec_file)->compound().at(
                jes_tag + "_L1L2L3Res_" + jec_algo);
        Logger::get("JetEnergyScaleData")
            ->debug("file: {}, function {}", jec_file,
                    (jes_tag + "_L1L2L3Res_" + jec_algo));
        auto JetEnergyScaleSF = [JES_evaluator](const float area,
                                                const float eta, const float pt,
                                                const float rho) {
            return JES_evaluator->evaluate({area, eta, pt, rho});
        };

        // lambda run with dataframe
        auto JetEnergyCorrectionLambda =
            [jes_tag,
             JetEnergyScaleSF](const ROOT::RVec<float> &pt_values,
                               const ROOT::RVec<float> &eta_values,
                               const ROOT::RVec<float> &area_values,
                               const ROOT::RVec<float> &rawFactor_values,
                               const float &rho_value) {
                ROOT::RVec<float> pt_values_corrected;
                for (int i = 0; i < pt_values.size(); i++) {
                    float corr_pt = pt_values.at(i);
                    if (jes_tag != "") {
                        // reapplying the JES correction
                        float raw_pt =
                            pt_values.at(i) * (1 - rawFactor_values.at(i));
                        float corr = JetEnergyScaleSF(area_values.at(i),
                                                      eta_values.at(i), raw_pt,
                                                      rho_value);
                        corr_pt = raw_pt * corr;
                        Logger::get("JetEnergyScaleData")
                            ->debug("reapplying JE scale for data: orig. jet "
                                    "pt {} to raw "
                                    "jet pt {} to recorr. jet pt {}",
                                    pt_values.at(i), raw_pt, corr_pt);
                    }
                    pt_values_corrected.push_back(corr_pt);
                    // if (pt_values_corrected.at(i)>15.0), this
                    // correction should be propagated to MET
                    // (requirement for type I corrections)
                }
                return pt_values_corrected;
            };
        auto df1 = df.Define(corrected_jet_pt, JetEnergyCorrectionLambda,
                             {jet_pt, jet_eta, jet_area, jet_rawFactor, rho});
        return df1;
    } else {
        auto df1 = df.Define(
            corrected_jet_pt,
            [](const ROOT::RVec<float> &pt_values) { return pt_values; },
            {jet_pt});
        return df1;
    }
}

/// Function to select jets passing a ID requirement, using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the rawID column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe
/// \param[in] idThreshold minimal ID value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutRawID(ROOT::RDF::RNode df, const std::string &quantity,
                          const std::string &maskname,
                          const float &idThreshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMin(idThreshold), {quantity});
    return df1;
}
/// Function to select jets failing a ID requirement, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the rawID column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe
/// \param[in] idThreshold maximal ID value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode AntiCutRawID(ROOT::RDF::RNode df, const std::string &quantity,
                              const std::string &maskname,
                              const float &idThreshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMax(idThreshold), {quantity});
    return df1;
}
} // end namespace jet
} // end namespace physicsobject

namespace quantities {
namespace jet {
/// Function to determine number of jets in a jet collection
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] jetcollection name of the vector that contains jet indices of the
/// jets belonging to the collection, its length constitutes the output quantity
///
/// \return a dataframe containing the number of jets in the jet collection
ROOT::RDF::RNode NumberOfJets(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &jetcollection) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &jetcollection) {
                         Logger::get("NumberOfJets")->debug("Counting jets");
                         Logger::get("NumberOfJets")
                             ->debug("NJets {}", jetcollection.size());
                         return (int)jetcollection.size();
                     },
                     {jetcollection});
}
/// Function to writeout the value of the btagger for a jet. The tag value is
/// identified via the discriminant of e.g. the DeepJet tagger from nanoAOD
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] btagcolumn name of the column that contains btag values of
/// the jets
/// \param[in] jetcollection name of the vector that contains jet indices of the
/// jets belonging to the collection, its length constitutes the output quantity
/// \param position The position in the jet collection vector, which is used to
/// store the index of the particle in the particle quantity vectors.
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode btagValue(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &btagcolumn,
                           const std::string &jetcollection,
                           const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &btagvalues,
                                const ROOT::RVec<int> &jetcollection) {
                         float btagValue = default_float;
                         try {
                             const int index = jetcollection.at(position);
                             btagValue = btagvalues.at(index);
                         } catch (const std::out_of_range &e) {
                         }
                         return btagValue;
                     },
                     {btagcolumn, jetcollection});
}
/// Function to writeout the hadron flavor for a jet.
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] jetcollection name of the vector that contains jet indices of the
/// jets belonging to the collection, its length constitutes the output quantity
/// \param position The position in the jet collection vector, which is used to
/// store the index of the particle in the particle quantity vectors.
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode flavor(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &flavorcolumn,
                        const std::string &jetcollection, const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &flavorvalues,
                                const ROOT::RVec<int> &jetcollection) {
                         int flavorValue = default_int;
                         const int index =
                             jetcollection.at(position, default_int);
                         flavorValue = flavorvalues.at(index, default_int);
                         return flavorValue;
                     },
                     {flavorcolumn, jetcollection});
}
} // end namespace jet
} // end namespace quantities
#endif /* GUARDJETS_H */
