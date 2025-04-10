#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "../include/SVFit/FastMTT.hxx"
#include "../include/SVFit/MeasuredTauLepton.hxx"
#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/vectoroperations.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
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


/// The namespace that is used to hold the functions for basic quantities that
/// are needed for every event
namespace quantities {
///mingxuan calc met MT2
// ROOT::RDF::RNode calc_MT2(ROOT::RDF::RNode df, const std::string &outputname,
//                         const std::string &mu1_p4, const std::string &mu2_p4, const std::string &met_p4) {
//     auto calc_MT2 = [](ROOT::Math::PtEtaPhiMVector &mu1_p4,
//                        ROOT::Math::PtEtaPhiMVector &mu2_p4,
//                        ROOT::Math::PtEtaPhiMVector &met_p4) {
//         float mVisA = mu1_p4.M();
//         float pxA = mu1_p4.Px();
//         float pyA = mu1_p4.Py();

//         float mVisB = mu2_p4.M();
//         float pxB = mu2_p4.Px();
//         float pyB = mu2_p4.Py();

//         float pxmiss = met_p4.Px();
//         float pymiss = met_p4.Py();

//         float chiA = 0;
//         float chiB = 0;

//         float desiredPrecisionOnMt2 = 0;

//         asymm_mt2_lester_bisect::disableCopyrightMessage();

//         float MT2 = asymm_mt2_lester_bisect::get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxmiss, pymiss, chiA, chiB, desiredPrecisionOnMt2);
//         if ( !std::isnan(MT2) && !std::isinf(MT2) ) {
//             return MT2;
//         } else {
//             return -10.0f;
//         }
//     };
//     return df.Define(outputname, calc_MT2, {mu1_p4, mu2_p4, met_p4});
// }
///mingxuan calc met Mct
ROOT::RDF::RNode calc_Mct(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &mu1_p4, const std::string &mu2_p4) {
    auto calc_met_Mct = [](ROOT::Math::PtEtaPhiMVector &mu1_p4,
                           ROOT::Math::PtEtaPhiMVector &mu2_p4) {
        float px1, px2, py1, py2, Et1, Et2;
        float Mct;
        px1 = mu1_p4.Px(), px2 = mu2_p4.Px(), py1 = mu1_p4.Py(), py2 = mu2_p4.Py();
        float muon_mass = 0.1056583755;
        Et1 = sqrt(px1*px1 + py1*py1 + muon_mass*muon_mass);
        Et2 = sqrt(px2*px2 + py2*py2 + muon_mass*muon_mass);
        Mct = sqrt((Et1+Et1)*(Et1+Et1) - ((px1-px2)*(px1-px2) + (py1-py2)*(py1-py2)));

        if ( !std::isnan(Mct) && !std::isinf(Mct) ) {
            return Mct;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calc_met_Mct, {mu1_p4, mu2_p4});
}
///mingxuan add A/B
ROOT::RDF::RNode calc_ratio(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &A, const std::string &B) {
    auto calculate_ratio = [](float A, float B) {
        if( B==0 ){
            return -10.0f;
        }else{
            return A/B;
        }
    };
    return df.Define(outputname, calculate_ratio, {A,B});
}
/// function to calculate the pz_nu
ROOT::RDF::RNode calculateNeutrinoPz(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &lep_p4, const std::string &met_p4) {
    auto calc_pz_nu = [](ROOT::Math::PtEtaPhiMVector &lep_p4,
                               ROOT::Math::PtEtaPhiMVector &met_p4) {
        const float mW = 80.379;
        float nu_px = met_p4.pt() * cos(met_p4.phi());
        float nu_py = met_p4.pt() * sin(met_p4.phi());
        //
        float a = (mW * mW) / 2 + lep_p4.Px() * nu_px + lep_p4.Py() * nu_py;
        float A = lep_p4.E() * lep_p4.E() - lep_p4.Pz() * lep_p4.Pz();
        float B = -2 * a * lep_p4.Pz();
        float C = lep_p4.E() * lep_p4.E() * (nu_px * nu_px + nu_py * nu_py) - a * a;
        // 
        float discriminant = B * B - 4 * A * C;
        float neutrino_pz1 = 0, neutrino_pz2 = 0;
        // 
        if (discriminant < 0) {
            // Complex solution: take real part
            neutrino_pz1 = -B / (2 * A);
            neutrino_pz2 = neutrino_pz1; // Identical since discriminant = 0
        } else {
            // Real solutions
            neutrino_pz1 = (-B + sqrt(discriminant)) / (2 * A);
            neutrino_pz2 = (-B - sqrt(discriminant)) / (2 * A);
        }
        if (abs(neutrino_pz1) < abs(neutrino_pz2)) {
            return neutrino_pz1;
        } else {
            return neutrino_pz2;
        }

    };
    return df.Define(outputname, calc_pz_nu, {lep_p4, met_p4});
}
///
/// function to pick dimuon pair from Higgs
///
///
ROOT::RDF::RNode calculate_kT(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p1, const std::string &p2) {
    auto calculate_kT = [](ROOT::Math::PtEtaPhiMVector &p1,
                               ROOT::Math::PtEtaPhiMVector &p2) {
        if (p1.pt() < 0.0 || p2.pt() < 0.0)
            return default_float;          
    float pt_1 = p1.pt();
    float pt_2 = p2.pt();
    float result =std::min(pt_1, pt_2) * ROOT::Math::VectorUtil::DeltaR(p1, p2);
    if ( !std::isnan(result) && !std::isinf(result) ) {
        return result;
    } else {
        return -10.0f;
    }
    };
    return df.Define(outputname, calculate_kT, {p1, p2});
}
ROOT::RDF::RNode calculate_antikT(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p1, const std::string &p2) {
    auto calculate_antikT = [](ROOT::Math::PtEtaPhiMVector &p1,
                               ROOT::Math::PtEtaPhiMVector &p2) {
        if (p1.pt() < 0.0 || p2.pt() < 0.0)
            return default_float;          
    float pt_1 = p1.pt();
    float pt_2 = p2.pt();
    float result =std::min(1/pt_1, 1/pt_2) * ROOT::Math::VectorUtil::DeltaR(p1, p2);
    if ( !std::isnan(result) && !std::isinf(result) ) {
        return result;
    } else {
        return -10.0f;
    }
    };
    return df.Define(outputname, calculate_antikT, {p1, p2});
}

///
/// funciton to calc mT contains MHT
ROOT::RDF::RNode mT_MHT(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met) {
    auto calculate_mt = [](ROOT::Math::PtEtaPhiMVector &particle_p4,
                           ROOT::Math::PtEtaPhiMVector &met) {
        float mt = (float)sqrt(2 * fabs(particle_p4.Pt()) * fabs(met.Pt()) *
                       (1. - cos(particle_p4.Phi() - met.Phi())));
        if ( !std::isnan(mt) && !std::isinf(mt) ) {
            return mt;
        } else {
            return default_float;
        }
    };
    return df.Define(outputname, calculate_mt, {particle_p4, met});
}
///
//// function to calc the delta Eta
ROOT::RDF::RNode deltaEta(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4) {
    auto calculate_deltaEta = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4) {
        return (float)fabs(p_1_p4.eta() - p_2_p4.eta());
    };
    return df.Define(outputname, calculate_deltaEta, {p_1_p4, p_2_p4});
}
///
/// Function to calculate the pt from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode pt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &inputvector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.pt(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode eta(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.eta(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode phi(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector) {
    return df.Define(outputname,
                     [](const ROOT::Math::PtEtaPhiMVector &p4) {
                         if (p4.pt() <
                             0.0) // negative pt is used to mark invalid LVs
                             return default_float;
                         return (float)p4.phi();
                     },
                     {inputvector});
}
/// Function to calculate the mass from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the mass value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode mass(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &inputvector) {
    return df.Define(outputname,
                     [](const ROOT::Math::PtEtaPhiMVector &p4) {
                         if (p4.pt() <
                             0.0) // negative pt is used to mark invalid LVs
                             return default_float;
                         return (float)p4.mass();
                     },
                     {inputvector});
}
/// Function to writeout the dxy impact parameter from a particle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the dxy value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param dxycolumn name of the column containing the dxy values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode dxy(ROOT::RDF::RNode df, const std::string &outputname,
                     const int &position, const std::string &pairname,
                     const std::string &dxycolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &dxy) {
            const int index = pair.at(position);
            return dxy.at(index, default_float);
        },
        {pairname, dxycolumn});
}
/// Function to writeout the dz impact parameter from a particle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the dz value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param dzcolumn name of the column containing the dz values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode dz(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &dzcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &dz) {
            const int index = pair.at(position);
            return dz.at(index, default_float);
        },
        {pairname, dzcolumn});
}
/// Function to writeout the charge of a particle. The particle is identified
/// via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the charge value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param chargecolumn name of the column containing the charge values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode charge(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &position, const std::string &pairname,
                        const std::string &chargecolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<int> &charge) {
            const int index = pair.at(position);
            return charge.at(index, default_int);
        },
        {pairname, chargecolumn});
}
/// Function to calculate the scalar sum of pts for given lorentz vectors and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode scalarPtSum(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &pt_1, const std::string &pt_2, const std::string &pt_3) {
    // build scalar sum of pts of 3 objects
    return df.Define(
        outputname,
        [](const float &pt_1,
           const float &pt_2, const float &pt_3) {
            if (pt_3 < 0.0 || pt_3 < 0.0 || pt_3 < 0.0)
                return default_float;
            auto const triple_lepton_pt = pt_1 + pt_2 + pt_3;
            return (float)triple_lepton_pt;
        },
        {pt_1, pt_2, pt_3});
}

ROOT::RDF::RNode scalarPtSum_met(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &pt_1, const std::string &pt_2, const std::string &pt_3 ,const std::string &pt_4 ) {
    // build scalar sum of pts of 3 objects
    return df.Define(
        outputname,
        [](const float &pt_1,
           const float &pt_2, const float &pt_3 , const float &pt_4) {
            if (pt_3 < 0.0 || pt_3 < 0.0 || pt_3 < 0.0 || pt_4 < 0.0)
                return default_float;
            auto const triple_lepton_pt_met = pt_1 + pt_2 + pt_3 + pt_4 ;
            return (float)triple_lepton_pt_met;
        },
        {pt_1, pt_2, pt_3, pt_4});
}

/**
 * @brief function used to calculate the deltaPhi between two lorentz vectors. $\phi_1$ is from the first lorentz vector and $\phi_2$ is from the second lorentz vector.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param p_1_p4 first lorentz vector
 * @param p_2_p4 second lorentz vector of
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4) {
    auto calculate_deltaPhi = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4) {
        return abs(ROOT::Math::VectorUtil::DeltaPhi(p_1_p4, p_2_p4));
    };
    return df.Define(outputname, calculate_deltaPhi, {p_1_p4, p_2_p4});
}
/**
 * @brief function used to calculate the deltaPhi between the lepton from a W and the visible Higgs decay products. $\phi_1$ is from the first lorentz vector and $\phi_2$ is from the second lorentz vector and \phi_3$ is from the third lorentz vector.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param p_1_p4 first lorentz vector
 * @param p_2_p4 second lorentz vector
 * @param p_3_p4 second lorentz vector
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaPhi_WH(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4, const std::string &p_3_p4) {
    auto calculate_deltaPhi = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4, ROOT::Math::PtEtaPhiMVector &p_3_p4) {
        auto const dileptonsystem = p_2_p4 + p_3_p4;
        return ROOT::Math::VectorUtil::DeltaPhi(p_1_p4, dileptonsystem);
    };
    return df.Define(outputname, calculate_deltaPhi, {p_1_p4, p_2_p4, p_3_p4});
}
/// Function to calculate the visible mass from a pair of lorentz vectors and
/// add it to the dataframe. The visible mass is calculated as the mass of the
/// lorentz vector of the dilepton system.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode m_vis(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::vector<std::string> &inputvectors) {
    // build visible mass from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.mass();
        },
        inputvectors);
}

/**
 * @brief Function used to calculate the FastMTT p4 from the given inputs. The
 * implementation is based on
 * https://github.com/SVfit/ClassicSVfit/tree/fastMTT_19_02_2019
 *
 * @param df The dataframe to add the quantity to
 * @param outputname name of the new column containing the lorentz vector value
 * @param pt_1 the name of the column containing the pt of the first particle
 * @param pt_2 the name of the column containing the pt of the second particle
 * @param eta_1  the name of the column containing the eta of the first particle
 * @param eta_2 the name of the column containing the eta of the second particle
 * @param phi_1 the name of the column containing the phi of the first particle
 * @param phi_2 the name of the column containing the phi of the second particle
 * @param mass_1 the name of the column containing the mass of the first
 * particle
 * @param mass_2 the name of the column containing the mass of the second
 * particle
 * @param met_pt the name of the column containing the met pt
 * @param met_phi the name of the column containing the met phi
 * @param met_cov_xx the name of the column containing the met covariance xx
 * @param met_cov_xy the name of the column containing the met covariance xy
 * @param met_cov_yy the name of the column containing the met covariance yy
 * @param decay_mode_1 the name of the column containing the decay mode of the
 * first particle
 * @param decay_mode_2 the name of the column containing the decay mode of the
 * second particle
 * @param finalstate the final state of the ditaudecay. Supported are "mt",
 * "et", "tt", "em"
 * @return ROOT::RDF::RNode
 */
// ROOT::RDF::RNode
// p4_fastmtt(ROOT::RDF::RNode df, const std::string &outputname,
//            const std::string &pt_1, const std::string &pt_2,
//            const std::string &eta_1, const std::string &eta_2,
//            const std::string &phi_1, const std::string &phi_2,
//            const std::string &mass_1, const std::string &mass_2,
//            const std::string &met_pt, const std::string &met_phi,
//            const std::string &met_cov_xx, const std::string &met_cov_xy,
//            const std::string &met_cov_yy, const std::string &decay_mode_1,
//            const std::string &decay_mode_2, const std::string &finalstate) {
//     auto calculate_fast_mtt =
//         [finalstate](const float &pt_1, const float &pt_2, const float &eta_1,
//                      const float &eta_2, const float &phi_1, const float &phi_2,
//                      const float &mass_1, const float &mass_2,
//                      const float &met_pt, const float &met_phi,
//                      const float &met_cov_xx, const float &met_cov_xy,
//                      const float &met_cov_yy, const int &decay_mode_1,
//                      const int &decay_mode_2) {
//             std::vector<fastmtt::MeasuredTauLepton> measuredTauLeptons;
//             TMatrixD covMET(2, 2);
//             covMET[0][0] = met_cov_xx;
//             covMET[1][0] = met_cov_xy;
//             covMET[0][1] = met_cov_xy;
//             covMET[1][1] = met_cov_yy;
//             // build the met lorentz vector
//             ROOT::Math::PtEtaPhiMVector met(met_pt, 0.0, met_phi, 0.0);
//             // set the decay modes according to the final state
//             auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
//             auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
//             int dm_1, dm_2;
//             if (finalstate == "mt") {
//                 dm_1 = -1;
//                 dm_2 = decay_mode_2;
//                 auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToMuDecay;
//                 auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
//             } else if (finalstate == "et") {
//                 dm_1 = -1;
//                 dm_2 = decay_mode_2;
//                 auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToElecDecay;
//                 auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
//             } else if (finalstate == "tt") {
//                 dm_1 = decay_mode_1;
//                 dm_2 = decay_mode_2;
//             } else if (finalstate == "em") {
//                 dm_1 = -1;
//                 dm_2 = -1;
//                 auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToElecDecay;
//                 auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToMuDecay;
//             } else {
//                 Logger::get("FastMTT")->error(
//                     "Final state {} not supported by FastMTT", finalstate);
//                 return (ROOT::Math::PtEtaPhiMVector)LorentzVector();
//             }
//             measuredTauLeptons.push_back(fastmtt::MeasuredTauLepton(
//                 decay_obj_1, pt_1, eta_1, phi_1, mass_1, dm_1));
//             measuredTauLeptons.push_back(fastmtt::MeasuredTauLepton(
//                 decay_obj_2, pt_2, eta_2, phi_2, mass_2, dm_2));
//             FastMTT FastMTTAlgo;
//             FastMTTAlgo.run(measuredTauLeptons, met.X(), met.Y(), covMET);
//             LorentzVector result = FastMTTAlgo.getBestP4();
//             // ROOT::Math::PtEtaPhiMVector result(_result.Pt(), _result.Eta(),
//             //                                    _result.Phi(), _result.M());
//             Logger::get("FastMTT")->debug("FastMTT result: {}", result.M());
//             return (ROOT::Math::PtEtaPhiMVector)result;
//         };
//     return df.Define(outputname, calculate_fast_mtt,
//                      {pt_1, pt_2, eta_1, eta_2, phi_1, phi_2, mass_1, mass_2,
//                       met_pt, met_phi, met_cov_xx, met_cov_xy, met_cov_yy,
//                       decay_mode_1, decay_mode_2});
// }
/// Function to calculate the visible pt from a pair of lorentz vectors and
/// add it to the dataframe. The visible pt is calculated as the pt of the
/// lorentz vector of the dilepton system.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pt_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.pt();
        },
        inputvectors);
}
/// Function to calculate the pt of the W from a the visible lepton fourvector, the met four vector and the neutrino four vector from the Higgs system and
/// add it to the dataframe.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pt_W(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2,
           const ROOT::Math::PtEtaPhiMVector &p4_3) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0 || p4_3.pt() < 0.0)
                return default_float;
            auto const w_p4 = p4_1 + p4_2 - p4_3;
            return (float)w_p4.Pt();
        },
        inputvectors);
}
/**
 * @brief Function to calculate the quantity `pZetaMissVis` from the two leptons
 in the event + the met vector. The variable is defined as
 * \f[
     D_\zeta = p_\zeta^\text{miss} - 0.85 p_\zeta^\text{vis}
    \qquad;
    p_\zeta^\text{miss} = \vec{p}_\text{T}^\text{miss} \cdot \hat{\zeta}
    \qquad;
    p_\zeta^\text{vis} = (\vec{p}_\text{T}^{p_1} + \vec{p}_\text{T}^{p_2}) \cdot
 \hat{\zeta} \f] where \f$\vec{p}_\text{T}^{p_{1,2}}\f$ corresponds to the
 transverse momentum vector of the first (second) lepton and \f$\hat{\zeta}\f$
 to the bisectional direction between the two leptons in the transverse plane.
 For more information check

 D. Jang, “Search for MSSM Higgs decaying to
 tau pairs in pp collision at √s=1.96 TeV at CDF”. PhD thesis, Rutgers
 University, 2006. FERMILAB-THESIS-2006-11.

 * @param df the input dataframe
 * @param p_1_p4 the lorentz vector of the first particle
 * @param p_2_p4 the lorentz vector of the second particle
 * @param met the lorentz vector of the met
 * @param outputname the name of the new column containing the pZetaMissVis
 value
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode pzetamissvis(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &p_1_p4,
                              const std::string &p_2_p4,
                              const std::string &met) {
    float alpha = 0.85;
    auto calculate_pzetamissvis = [alpha](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                                          ROOT::Math::PtEtaPhiMVector &p_2_p4,
                                          ROOT::Math::PtEtaPhiMVector &met) {
        auto met_3dvec = met.Vect();
        met_3dvec.SetZ(0.0);
        // calculate zeta for the delepton system
        auto p1_norm = p_1_p4.Vect().Unit();
        auto p2_norm = p_2_p4.Vect().Unit();
        p1_norm.SetZ(0.0);
        p2_norm.SetZ(0.0);
        p1_norm = p1_norm.Unit();
        p2_norm = p2_norm.Unit();
        auto zeta = (p1_norm + p2_norm).Unit();

        auto dileptonsystem = p_1_p4.Vect() + p_2_p4.Vect();
        dileptonsystem.SetZ(0);
        auto pzetaVis = dileptonsystem.Dot(zeta);
        return met_3dvec.Dot(zeta) - (alpha * pzetaVis);
    };
    return df.Define(outputname, calculate_pzetamissvis, {p_1_p4, p_2_p4, met});
}
/**
 * @brief function used to calculate mTdileptonMET, which is the transverse mass
 * of the di-lepton system. The transverse mass is calculated using the
 * vectoroperations::calculateMT function.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the mTdileptonMET value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode mTdileptonMET(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &p_1_p4,
                               const std::string &p_2_p4,
                               const std::string &met) {
    auto calculate_mTdileptonMET = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                                      ROOT::Math::PtEtaPhiMVector &p_2_p4,
                                      ROOT::Math::PtEtaPhiMVector &met) {
        ROOT::Math::PtEtaPhiMVector dilepton = p_1_p4 + p_2_p4;
        return vectoroperations::calculateMT(dilepton, met);
    };
    return df.Define(outputname, calculate_mTdileptonMET,
                     {p_1_p4, p_2_p4, met});
}

/**
 * @brief function used to calculate the deltaR between two lorentz vectors. It
 is defined as \f[ \Delta R = \sqrt{(\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2}
 \f$ where $\eta_1$ and $\phi_1$ are from the first lorentz vector and $\eta_2$
 and $\phi_2$ are from the second lorentz vector.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param p_1_p4 first lorentz vector
 * @param p_2_p4 second lorentz vector of
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaR(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4) {
    auto calculate_deltaR = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4) {
        return (float)ROOT::Math::VectorUtil::DeltaR(p_1_p4, p_2_p4);
    };
    return df.Define(outputname, calculate_deltaR, {p_1_p4, p_2_p4});
}
/**
 * @brief function used to calculate the transverse mass of a particle. The
 * transverse mass is calculated using the vectoroperations::calculateMT
 * function.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the mT value
 * @param particle_p4 lorentz vector of the particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode mT(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met) {
    auto calculate_mt = [](ROOT::Math::PtEtaPhiMVector &particle_p4,
                           ROOT::Math::PtEtaPhiMVector &met) {
        return vectoroperations::calculateMT(particle_p4, met);
    };
    return df.Define(outputname, calculate_mt, {particle_p4, met});
}

/**
 * @brief function used to calculate the pt of the dilepton + met system.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the pt_tt value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode pt_tt(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &p_1_p4, const std::string &p_2_p4,
                       const std::string &met) {
    auto calculate_pt_tt = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                              ROOT::Math::PtEtaPhiMVector &p_2_p4,
                              ROOT::Math::PtEtaPhiMVector &met) {
        auto dileptonmet = p_1_p4 + p_2_p4 + met;
        return (float)dileptonmet.Pt();
    };
    return df.Define(outputname, calculate_pt_tt, {p_1_p4, p_2_p4, met});
}

/**
 * @brief function used to calculate the pt of the dilepton + two leading jets +
 * met system. If the number of jets is less than 2, the quantity is set to 10
 * instead.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the pt_ttjj value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param jet_1_p4 lorentz vector of the first jet
 * @param jet_2_p4 lorentz vector of the second jet
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode pt_ttjj(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::string &p_1_p4, const std::string &p_2_p4,
                         const std::string &jet_1_p4,
                         const std::string &jet_2_p4, const std::string &met) {
    auto calculate_pt_ttjj = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                                ROOT::Math::PtEtaPhiMVector &p_2_p4,
                                ROOT::Math::PtEtaPhiMVector &jet_1_p4,
                                ROOT::Math::PtEtaPhiMVector &jet_2_p4,
                                ROOT::Math::PtEtaPhiMVector &met) {
        if (jet_1_p4.pt() < 0.0 || jet_2_p4.pt() < 0.0)
            return default_float;
        auto jetlepmet = p_1_p4 + p_2_p4 + met + jet_1_p4 + jet_2_p4;
        return (float)jetlepmet.Pt();
    };
    return df.Define(outputname, calculate_pt_ttjj,
                     {p_1_p4, p_2_p4, jet_1_p4, jet_2_p4, met});
}

/**
 * @brief function used to calculate the pt two leading jets
 If the number of jets is less than 2, the quantity is set to 10
 * instead.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the pt_dijet value
 * @param jet_1_p4 lorentz vector of the first jet
 * @param jet_2_p4 lorentz vector of the second jet
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode pt_dijet(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &jet_1_p4,
                          const std::string &jet_2_p4) {
    auto calculate_pt_dijet = [](ROOT::Math::PtEtaPhiMVector &jet_1_p4,
                                 ROOT::Math::PtEtaPhiMVector &jet_2_p4) {
        if (jet_1_p4.pt() < 0.0 || jet_2_p4.pt() < 0.0)
            return default_float;
        auto dijetsystem = jet_1_p4 + jet_2_p4;
        return (float)dijetsystem.Pt();
    };
    return df.Define(outputname, calculate_pt_dijet, {jet_1_p4, jet_2_p4});
}

/**
 * @brief function used to check the hemisphere of the two leading jets, if both
 jets are in the same hemisphere, the quantity is set to 1, otherwise it is set
 to 0. If the number of jets is less than 2, the quantity is set to 10
 * instead.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the jet hemisphere value
 * @param jet_1_p4 lorentz vector of the first jet
 * @param jet_2_p4 lorentz vector of the second jet
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode jet_hemisphere(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &jet_1_p4,
                                const std::string &jet_2_p4) {
    auto calculate_jet_hemisphere = [](ROOT::Math::PtEtaPhiMVector &jet_1_p4,
                                       ROOT::Math::PtEtaPhiMVector &jet_2_p4) {
        if (jet_1_p4.pt() < 0.0 || jet_2_p4.pt() < 0.0)
            return default_int;
        return (int)(jet_1_p4.Eta() * jet_2_p4.Eta() > 0);
    };
    return df.Define(outputname, calculate_jet_hemisphere,
                     {jet_1_p4, jet_2_p4});
}

/**
 * @brief function used to calculate `mt_tot`. It is defined as
 \f[
    m_{T}^{tot} = \sqrt{m_{T}^2(p_{1},E_{T}^{miss}) +
 m_{T}^2(p_{2},E_{T}^{miss}) + m_{T}^2(p_{1},p_2) } \f] where \f$ m_{T}^2 \f$ is
 the transverse mass, \f$ p_{1} \f$ and \f$ p_{2} \f$ are the lepton lorentz
 vectors and \f$ E_{T}^{miss} \f$ is the missing energy.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the mt_tot value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode mt_tot(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met) {
    auto calculate_mt_tot = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4,
                               ROOT::Math::PtEtaPhiMVector &met) {
        const float mt_1 = vectoroperations::calculateMT(p_1_p4, met);
        const float mt_2 = vectoroperations::calculateMT(p_2_p4, met);
        const float mt_mix = vectoroperations::calculateMT(p_1_p4, p_2_p4);
        return (float)sqrt(mt_1 * mt_1 + mt_2 * mt_2 + mt_mix * mt_mix);
    };
    return df.Define(outputname, calculate_mt_tot, {p_1_p4, p_2_p4, met});
}

/// Function to writeout the isolation of a particle. The particle is
/// identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the isolation value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param isolationcolumn name of the column containing the isolation
/// values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode isolation(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &isolationcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<float> &isolation) {
                         const int index = pair.at(position);
                         return isolation.at(index, default_float);
                     },
                     {pairname, isolationcolumn});
}
/// Function to writeout the PDGID from a genparticle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the PDGID
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param pdgidcolumn name of the column containing the pdgID
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pdgid(ROOT::RDF::RNode df, const std::string &outputname,
                       const int &position, const std::string &pairname,
                       const std::string &pdgidcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<int> &pdgid) {
            const int index = pair.at(position);
            return pdgid.at(index, default_pdgid);
        },
        {pairname, pdgidcolumn});
}

// vhmm extend to number of objects, not only leptons
/// Function to determine number of good objects
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] goodleptons name of the vector that contains a lepton mask of
/// good leptons, its length of non-zero values constitutes the output quantity
///
/// \return a dataframe containing the number of good leptons in an event
ROOT::RDF::RNode NumberOfGoodObjects(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const std::string &goodmask) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &goodmask) {
                         return (int)ROOT::VecOps::Nonzero(goodmask).size();
                     },
                     {goodmask});
}

/// Function to determine number of good leptons
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] goodleptons name of the vector that contains a lepton mask of
/// good leptons, its length of non-zero values constitutes the output quantity
///
/// \return a dataframe containing the number of good leptons in an event
ROOT::RDF::RNode NumberOfGoodLeptons(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const std::string &goodleptons) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &goodleptons) {
                         return (int)ROOT::VecOps::Nonzero(goodleptons).size();
                     },
                     {goodleptons});
}
/// namespace for tau specific quantities
namespace tau {
/// Function to writeout the decaymode of a tau. The particle is identified
/// via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the decaymode value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param decaymodecolumn name of the column containing the decaymode
/// values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode decaymode(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &decaymodecolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &decaymode) {
                         const int index = pair.at(position);
                         return decaymode.at(index, default_int);
                     },
                     {pairname, decaymodecolumn});
}
/// Function to writeout the genmatch of a tau. The particle is identified
/// via the index stored in the pair vector Genmatch values are defined as
/// \code
///   1 = prompt electron,
///   2 = prompt muon,
///   3 = tau->e decay,
///   4 = tau->mu decay,
///   5 = hadronic tau decay,
///   0 = unknown or unmatched
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the genmatch value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param genmatchcolumn name of the column containing the genmatch values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode genmatch(ROOT::RDF::RNode df, const std::string &outputname,
                          const int &position, const std::string &pairname,
                          const std::string &genmatchcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<UChar_t> &genmatch) {
                         const int index = pair.at(position);
                         return genmatch.at(index, default_uchar);
                     },
                     {pairname, genmatchcolumn});
}
/// Function to writeout the pt of the reco jet associated with the given
/// tau.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association
/// between the tau and the reco jet \param jetpt_column name of the column
/// containing the recojet pt values
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode matching_jet_pt(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const int &position,
                                 const std::string &pairname,
                                 const std::string &taujet_index,
                                 const std::string &jetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<float> &jetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         return jetpt.at(jetindex, default_float);
                     },
                     {pairname, taujet_index, jetpt_column});
}
/// Function to writeout the pt of the gen jet associated with the reco jet,
/// which is associated with the given tau. \code
///  Tau --> recoJet --> GenJet
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association
/// between the tau and the reco jet \param genjet_index name of the column
/// containing the association between the reco jet and the gen jet \param
/// genjetpt_column name of the column containing the genJet pt values
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode matching_genjet_pt(
    ROOT::RDF::RNode df, const std::string &outputname, const int &position,
    const std::string &pairname, const std::string &taujet_index,
    const std::string &genjet_index, const std::string &genjetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<int> &genjets,
                                const ROOT::RVec<float> &genjetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         const int genjetindex = genjets.at(jetindex, -1);
                         return genjetpt.at(genjetindex, default_float);
                     },
                     {pairname, taujet_index, genjet_index, genjetpt_column});
}
/// Function to writeout a flag if a tau passes a specific tau id cut. The
/// particle is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the flag
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param nameID name of the ID column in the NanoAOD
/// \param idxID bitvalue of the WP the has to be passed
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode TauIDFlag(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &nameID, const int &idxID) {
    return df.Define(
        outputname,
        [position, idxID](const ROOT::RVec<int> &pair,
                          const ROOT::RVec<UChar_t> &IDs) {
            Logger::get("tauIDFlag")
                ->debug(
                    "position tau in pair {}, pair {}, id bit {}, vsjet ids {}",
                    position, pair, idxID, IDs);
            const int index = pair.at(position);
            const int ID = IDs.at(index, default_int);
            if (ID != default_int)
                return std::min(1, int(ID & 1 << (idxID - 1)));
            else
                return int(ID);
        },
        {pairname, nameID});
}
} // end namespace tau
/// namespace for muon specific quantities
namespace muon {
/**
 * @brief Function to writeout the id of a muon.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param idcolumn name of the column containing the muon id values
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<bool> &id) {
            const int index = pair.at(position, -1);
            return id.at(index, false);
        },
        {pairname, idcolumn});
}
/**
 * @brief Function to writeout the is global flag of a muon.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param globalflagcolumn name of the column containing the muon is global flag
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode is_global(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &globalflagcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<bool> &globalflag) {
                         const int index = pair.at(position, -1);
                         return globalflag.at(index, false);
                     },
                     {pairname, globalflagcolumn});
}
/**
 * @brief Function to writeout the is tracker flag of a muon.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param trackerflagcolumn name of the column containing the muon is global flag
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode is_tracker(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &trackerflagcolumn) {
    Logger::get("muonIsTrackerflag")
                ->debug(
                    "is tracker pos {}", position);
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<bool> &trackerflag) {
                         const int index = pair.at(position, -1);
                         return trackerflag.at(index, false);
                     },
                     {pairname, trackerflagcolumn});
}
} // end namespace muon
namespace electron {
/**
 * @brief Function to writeout the id of a electron.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param idcolumn name of the column containing the muon id values
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn) {
    Logger::get("electronIDflag")
                ->debug(
                    "ele ID position {}", position);
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<bool> &id) {
            const int index = pair.at(position, -1);
            return id.at(index, false);
        },
        {pairname, idcolumn});
}
} // end namespace electron

namespace boostedbbtt{
    ROOT::RDF::RNode dR_tau_fatjet0(
        ROOT::RDF::RNode df, 
        const std::string &dR_tau_fatjet,
        const std::string &p0_fatjet, 
        const std::string &taupt, 
        const std::string &taueta, 
        const std::string &tauphi,
        const std::string &taumass
    ) {
        auto df1 = df.Define(
            dR_tau_fatjet,
            [](
                const ROOT::Math::PtEtaPhiMVector &p0,
                const ROOT::RVec<float> &pt, 
                const ROOT::RVec<float> &eta, 
                const ROOT::RVec<float> &phi, 
                const ROOT::RVec<float> &mass 
            ) {
                ROOT::RVec<float> deltaRs;

                for (size_t i = 0; i < pt.size(); ++i) {
                    ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], mass[i]);
                    float deltaR_value = (float)ROOT::Math::VectorUtil::DeltaR(p0, p4);
                    deltaRs.push_back(deltaR_value);
                    }
                    return deltaRs;
                },
                {p0_fatjet ,taupt, taueta, tauphi, taumass}
        );

        return df1;
    }

    void DecomposeMomentum(
        const ROOT::Math::PtEtaPhiMVector& p, 
        const ROOT::Math::PtEtaPhiMVector& a, 
        const ROOT::Math::PtEtaPhiMVector& b,
        ROOT::Math::PtEtaPhiMVector& pa,
        ROOT::Math::PtEtaPhiMVector& pb
    ) {
        //
        typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> CartesianVector;

        //
        CartesianVector p_xyz(p.Px(), p.Py(), p.Pz(), p.E());
        CartesianVector a_xyz(a.Px(), a.Py(), a.Pz(), a.E());
        CartesianVector b_xyz(b.Px(), b.Py(), b.Pz(), b.E());

        //
        double a_magnitude = a_xyz.R();
        double b_magnitude = b_xyz.R();

        // 
        CartesianVector a_unit = a_xyz / a_magnitude;
        CartesianVector b_unit = b_xyz / b_magnitude;

        //
        double p_dot_a = p_xyz.Dot(a_unit);
        double p_dot_b = p_xyz.Dot(b_unit);

        CartesianVector p_a = p_dot_a * a_unit;
        CartesianVector p_b = p_dot_b * b_unit;

        //
        pa = ROOT::Math::PtEtaPhiMVector(p_a.Pt(), p_a.Eta(), p_a.Phi(), p_a.M());
        pb = ROOT::Math::PtEtaPhiMVector(p_b.Pt(), p_b.Eta(), p_b.Phi(), p_b.M());
    }

    ROOT::RDF::RNode Mass_CA(
        ROOT::RDF::RNode df,
        const std::string &p0_fatjet,
        const std::string &p1_fatjet,
        const std::string &met_pt,
        const std::string &met_phi,
        const std::string &mass_caca
    ) 
    {
        auto df1 = df.Define(
            mass_caca,
            [](
                const ROOT::Math::PtEtaPhiMVector &p0,
                const ROOT::Math::PtEtaPhiMVector &p1,
                const float &pt_met,
                const float &phi_met) {

                    ROOT::Math::PtEtaPhiMVector p_MET(pt_met, 0, phi_met, pt_met);
                    ROOT::Math::PtEtaPhiMVector p_Nu0, p_Nu1;
                    DecomposeMomentum(p_MET, p0, p1, p_Nu0, p_Nu1);

                    float Nu0Tau_pt = p_Nu0.Pt();
                    float Nu1Tau_pt = p_Nu1.Pt();

                    float x0_METtoNu_tau = p0.Pt() / (p0.Pt() + Nu0Tau_pt);
                    float x1_METtoNu_tau = p1.Pt() / (p1.Pt() + Nu1Tau_pt);

                    float x0x1 = x0_METtoNu_tau * x1_METtoNu_tau;

                    return x0x1;
                    },
                    {p0_fatjet, p1_fatjet, met_pt, met_phi}
        );
        return df1;
    }

    ROOT::RDF::RNode CA_ttMAss(
        ROOT::RDF::RNode df,
        const std::string &x0x1,
        const std::string &FatJet_mass,
        const std::string &result
    ) 
    {
        auto df1 = df.Define(
            result,
            [](
                const float &xx,
                const float &mass) {
                    float Mass = mass/std::sqrt(xx);
                    return Mass;
                    },
                    {x0x1, FatJet_mass}
        );
        return df1;
    }

    ROOT::RDF::RNode CA_ttMAss_fatjet(
        ROOT::RDF::RNode df,
        const std::string &p0_fatjet,
        const std::string &met_pt,
        const std::string &met_phi,
        const std::string &result
    ) 
    {
        auto df1 = df.Define(
            result,
            [](
                const ROOT::Math::PtEtaPhiMVector &p0,
                const float &pt_met,
                const float &phi_met) {
                    ROOT::Math::PtEtaPhiMVector p_MET(pt_met, 0, phi_met, pt_met);
                    ROOT::Math::PtEtaPhiMVector p_Nu0, p_Nu1;
                    DecomposeMomentum(p_MET, p0, p0, p_Nu0, p_Nu1);
                    float mm = (p0+p_Nu0).M();
                    return mm;
                    },
                    {p0_fatjet, met_pt, met_phi}
        );
        return df1;
    }

    ROOT::RDF::RNode Score_9X(
        ROOT::RDF::RNode df,
        const std::string &Score,
        const std::string &result
    ) 
    {
        auto df1 = df.Define(
            result,
            [](const float &score) 
            {
                    float mm = -std::log10(1 - score + 1e-18);
                    return mm;
                    },
                    {Score}
        );
        return df1;
    }

    ROOT::RDF::RNode Score_original(
        ROOT::RDF::RNode df,
        const std::string &Score,
        const std::string &QCD,
        const std::string &result
    ) 
    {
        auto df1 = df.Define(
            result,
            [](
                const float &score,
                const float &qcd) 
            {
                    float mm = qcd * score / (1-score);
                    return mm;
                    },
                    {Score, QCD}
        );
        return df1;
    }

}

} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */
