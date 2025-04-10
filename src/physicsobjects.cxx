#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

#include "../include/defaults.hxx"
#include "../include/RoccoR.hxx"
#include "../include/basefunctions.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "correction.h"
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
/// Namespace containing function to apply cuts on physics objects. The
/// cut results are typically stored within a mask, which is represented by
/// an `ROOT::RVec<int>`.
///    \code
///    In the mask
///    1 --> cut is passed by the object
///    0 --> cut is not passed by the object
///    \endcode
/// multiple cuts can be combined by multiplying masks using
/// physicsobject::CombineMasks.
namespace physicsobject {
///mingxuan calc 4l cosThStar from HZZ
ROOT::RDF::RNode calc_4l_cosThStar_hzz(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &v_p4_, const std::string &h_p4_, const std::string &mu_p4_) {
    auto calc_cosThStar = [](ROOT::Math::PtEtaPhiMVector &v_p4_,
                             ROOT::Math::PtEtaPhiMVector &h_p4_,
                             ROOT::Math::PtEtaPhiMVector &mu_p4_) {
        TLorentzVector h_p4;
        TLorentzVector v_p4;
        TLorentzVector mu_p4;
        h_p4.SetPtEtaPhiM(h_p4_.Pt(), h_p4_.Eta(), h_p4_.Phi(), h_p4_.M());
        v_p4.SetPtEtaPhiM(v_p4_.Pt(), v_p4_.Eta(), v_p4_.Phi(), v_p4_.M());
        mu_p4.SetPtEtaPhiM(mu_p4_.Pt(), mu_p4_.Eta(), mu_p4_.Phi(), mu_p4_.M());

        TLorentzVector mom_v_p4;
        mom_v_p4 = v_p4 + h_p4;
        TVector3 boost_vec;
        boost_vec = h_p4.BoostVector();
        boost_vec = -boost_vec;
        mu_p4.Boost(boost_vec);
        TVector3 mu_vec = mu_p4.Vect();
        TVector3 mom_v_vec = mom_v_p4.Vect();
        float cosThStar_4l_hzz = cos(mu_vec.Angle(mom_v_vec));

        if ( !std::isnan(cosThStar_4l_hzz) && !std::isinf(cosThStar_4l_hzz) ) {
            return cosThStar_4l_hzz;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calc_cosThStar, {v_p4_, h_p4_, mu_p4_});
}
///mingxuan calc 4l costheta1 from hzz
ROOT::RDF::RNode calc_4l_cosTh1_hzz(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &v_p4_, const std::string &h_p4_) {
    auto calc_cosTh1 = [](ROOT::Math::PtEtaPhiMVector &v_p4_,
                          ROOT::Math::PtEtaPhiMVector &h_p4_) {
        TLorentzVector h_p4;
        TLorentzVector v_p4;
        h_p4.SetPtEtaPhiM(h_p4_.Pt(), h_p4_.Eta(), h_p4_.Phi(), h_p4_.M());
        v_p4.SetPtEtaPhiM(v_p4_.Pt(), v_p4_.Eta(), v_p4_.Phi(), v_p4_.M());
        
        TLorentzVector mom_v_p4;
        mom_v_p4 = v_p4 + h_p4;
        TVector3 quark(0, 0, 1);
        TVector3 mom_v_vec = mom_v_p4.Vect();
        float cosTh1_hzz = cos(quark.Angle(mom_v_vec));

        if ( !std::isnan(cosTh1_hzz) && !std::isinf(cosTh1_hzz) ) {
            return cosTh1_hzz;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calc_cosTh1, {v_p4_, h_p4_});
}
///mingxuan calc 4l cosphi from hzz
ROOT::RDF::RNode calc_4l_cosphi_hzz(ROOT::RDF::RNode df, const std::string &outputname,
                                    const std::string &lep_v4_, const std::string &v_p4_, const std::string &h_p4_) {
    auto calc_cosphi = [](ROOT::Math::PtEtaPhiMVector &lep_p4_,
                          ROOT::Math::PtEtaPhiMVector &v_p4_,
                          ROOT::Math::PtEtaPhiMVector &h_p4_){
        TLorentzVector h_p4;
        TLorentzVector v_p4;
        TLorentzVector lep_p4;
        h_p4.SetPtEtaPhiM(h_p4_.Pt(), h_p4_.Eta(), h_p4_.Phi(), h_p4_.M());
        v_p4.SetPtEtaPhiM(v_p4_.Pt(), v_p4_.Eta(), v_p4_.Phi(), v_p4_.M());
        lep_p4.SetPtEtaPhiM(lep_p4_.Pt(), lep_p4_.Eta(), lep_p4_.Phi(), lep_p4_.M());
        
        TLorentzVector mom_v_p4;
        mom_v_p4 = v_p4 + h_p4;
        TVector3 boost_vec;
        boost_vec = v_p4.BoostVector();
        lep_p4.Boost(-boost_vec);
        TVector3 flat1;
        TVector3 flat2;
        TVector3 quark(0, 0, 1);
        flat1 = (lep_p4.Vect()).Cross(v_p4.Vect()), flat2 = (mom_v_p4.Vect()).Cross(quark);
        float cosphi = cos(flat1.Angle(flat2));

        if ( !std::isnan(cosphi) && !std::isinf(cosphi) ) {
            return cosphi;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calc_cosphi, {lep_v4_, v_p4_, h_p4_});
}
///mingxuan calc 4l cosphi1 from hzz
ROOT::RDF::RNode calc_4l_cosphi1_hzz(ROOT::RDF::RNode df, const std::string &outputname,
                                     const std::string &mu_p4_, const std::string &h_p4_, const std::string &v_p4_) {
    auto calc_cosphi1 = [](ROOT::Math::PtEtaPhiMVector &mu_p4_, 
                           ROOT::Math::PtEtaPhiMVector &h_p4_,
                           ROOT::Math::PtEtaPhiMVector &v_p4_) {
        TLorentzVector h_p4;
        TLorentzVector v_p4;
        TLorentzVector mu_p4;
        h_p4.SetPtEtaPhiM(h_p4_.Pt(), h_p4_.Eta(), h_p4_.Phi(), h_p4_.M());
        v_p4.SetPtEtaPhiM(v_p4_.Pt(), v_p4_.Eta(), v_p4_.Phi(), v_p4_.M());
        mu_p4.SetPtEtaPhiM(mu_p4_.Pt(), mu_p4_.Eta(), mu_p4_.Phi(), mu_p4_.M());
        
        TLorentzVector mom_v_p4;
        mom_v_p4 = v_p4 + h_p4;
        TVector3 boost_vec;
        boost_vec = h_p4.BoostVector();
        mu_p4.Boost(-boost_vec);
        TVector3 flat1;
        TVector3 flat2;
        TVector3 quark(0, 0, 1);
        flat1 = (mu_p4.Vect()).Cross(mom_v_p4.Vect()), flat2 = (mom_v_p4.Vect()).Cross(quark);
        float cosphi1 = cos(flat1.Angle(flat2));

        if ( !std::isnan(cosphi1) && !std::isinf(cosphi1) ) {
            return cosphi1;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calc_cosphi1, {mu_p4_, h_p4_, v_p4_});
}

/// write by botao
///
/// function to calculate the W p4
ROOT::RDF::RNode build_Wp4(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &lep_p4, const std::string &nu_p4) {
    auto calc_Wp4 = [](ROOT::Math::PtEtaPhiMVector &lep_p4, ROOT::Math::PtEtaPhiMVector &nu_p4) {
        ROOT::Math::PtEtaPhiMVector W_p4;
        W_p4 = lep_p4 + nu_p4;
        return W_p4;
    };
    auto df1 = 
        df.Define(outputname, calc_Wp4, {lep_p4, nu_p4});
    return df1;
}
///
/// function to calculate the nu p4
ROOT::RDF::RNode build_nup4(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &met_p4, const std::string &pz_nu) {
    auto calc_nup4 = [](ROOT::Math::PtEtaPhiMVector &met_p4, const float &pz_nu) {
        TLorentzVector nu_tmp_p4;
        ROOT::Math::PtEtaPhiMVector nu_p4;
        // ROOT::Math::PxPyPzEVector nu_p4;
        float nu_px = met_p4.pt() * cos(met_p4.phi());
        float nu_py = met_p4.pt() * sin(met_p4.phi());
        float nu_e = sqrt(nu_px * nu_px + nu_py * nu_py + pz_nu * pz_nu);
        nu_tmp_p4.SetPxPyPzE(nu_px, nu_py, pz_nu, nu_e);
        nu_p4 = ROOT::Math::PtEtaPhiMVector(nu_tmp_p4.Pt(), nu_tmp_p4.Eta(), nu_tmp_p4.Phi(), nu_tmp_p4.M());
        // nu_p4 = ROOT::Math::PxPyPzEVector(nu_px, nu_py, pz_nu, nu_e);

        return nu_p4;
    };
    auto df1 = 
        df.Define(outputname, calc_nup4, {met_p4, pz_nu});
    return df1;
}
///
/// function to pick dimuon pair from Higgs
ROOT::RDF::RNode Muon_var(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_vars,
                                 const std::string &dimuons_index,
                                 const int ranking) {
    auto muon_var = [ranking](const ROOT::RVec<float> &muon_vars,
                                 const ROOT::RVec<int> &dimuons_index) {
                                 if ( dimuons_index.at(ranking) == -1 ) {
                                    return default_float;
                                 } else {
                                    // ranking should be 0 or 1, stands for leading and sub-leading muon
                                    return (float)muon_vars.at(dimuons_index[ranking]);
                                 }
                             };
    auto df1 = 
        df.Define(outputname, muon_var, {muon_vars, dimuons_index});
    return df1;
}
///
/// funciton to calc pt_W using extra_lep_p4 and met_p4
ROOT::RDF::RNode pt_W(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met) {
    auto pt_w = [](ROOT::Math::PtEtaPhiMVector &particle_p4,
                           ROOT::Math::PtEtaPhiMVector &met) {
        if ( particle_p4.Pt() < 0 || met.Pt() < 0 )
            return default_float;
        auto const W_system = particle_p4 + met;
        return (float)W_system.Pt();
    };
    return df.Define(outputname, pt_w, {particle_p4, met});
}
///
/// funciton to calc phi_W using extra_lep_p4 and met_p4
ROOT::RDF::RNode phi_W(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met) {
    auto phi_w = [](ROOT::Math::PtEtaPhiMVector &particle_p4,
                           ROOT::Math::PtEtaPhiMVector &met) {
        if ( particle_p4.Pt() < 0 || met.Pt() < 0 )
            return default_float;
        auto const W_system = particle_p4 + met;
        return (float)W_system.Phi();
    };
    return df.Define(outputname, phi_w, {particle_p4, met});
}
///

/// function to set events that jet_pt=-999 -> flag=0
ROOT::RDF::RNode PassJetVetoFlag(ROOT::RDF::RNode df, const std::string &jet_pts, const std::string &outputname) {
    auto PassJetVetoFlag = [](const ROOT::RVec<float> &jet_pts) {
        /// loop the jet pt to find if any == -999
        for (unsigned int j = 0; j < (int)jet_pts.size(); ++j ) {
            if (jet_pts.at(j) == -999.0) {
                return 0;
            }
        }
        return 1;
    };
    auto df1 = 
        df.Define(outputname, PassJetVetoFlag, {jet_pts});
    return df1;
}
///

/// 
ROOT::RDF::RNode flagNumObject(ROOT::RDF::RNode df, const std::string &flagname,
                               const std::string &object_number,
                               const int flag_threshold, const std::string relation) {
    // std::cout << "Check point" << std::endl;
    auto df1 = 
        df.Define( flagname, object_number + " " + relation + " " + std::to_string(flag_threshold) );
    return df1;
}

/// function to select the smallest mass of dilepton pair
ROOT::RDF::RNode M_dileptonMass(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &particle_charges,
                                 const std::string &goodmuons_index) {
    auto mass_calculation = [](const ROOT::RVec<float> &particle_pts,
                               const ROOT::RVec<float> &particle_etas,
                               const ROOT::RVec<float> &particle_phis,
                               const ROOT::RVec<float> &particle_masses,
                               const ROOT::RVec<int> &particle_charges,
                               const ROOT::RVec<int> &goodmuons_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    try {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodmuons_index[k]), 
                                                                         particle_etas.at(goodmuons_index[k]),
                                                                         particle_phis.at(goodmuons_index[k]),
                                                                         particle_masses.at(goodmuons_index[k])));
                                    } catch (const std::out_of_range &e) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_1;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_2;
                                 p4_1 = p4;
                                 p4_2 = p4;
                                 std::vector<float> masses;
                                 for (unsigned int i = 0; i < p4_1.size(); ++i) {
                                     for (unsigned int j = i + 1; j < p4_2.size(); ++j) {
                                         if (p4_1[i].pt() < 0.0 || p4_2[j].pt() < 0.0) {
                                             continue;
                                         }
                                         if ( particle_charges[goodmuons_index[i]] + particle_charges[goodmuons_index[j]] != 0 ) {
                                             continue;
                                         }
                                         auto dileptonsystem = p4_1[i] + p4_2[j];
                                         masses.push_back((float)dileptonsystem.mass());
                                     }
                                 }
                                 std::sort(masses.begin(), masses.end());
                                 if ( masses.size() <= 0 ) {
                                    return 999.0f;
                                 }
                                 else {
                                    return masses[0];
                                 }
                             };
    // std::vector<std::string> column_names = {particle_pts, particle_etas, particle_phis, particle_masses, muon_size};
    auto df1 = 
        df.Define(outputname, mass_calculation, {particle_pts, particle_etas, particle_phis, particle_masses, particle_charges, goodmuons_index});
    return df1;
}

/// function to veto ECal Gap
/// \param[in] etaBoundary The upper boundary of the eta, such as 2.5
///
ROOT::RDF::RNode ECalGapVeto(ROOT::RDF::RNode df, const std::string &etaColumnName,
                              const std::string &maskname,
                              const float &etaBoundary, const float &lowerThresholdBarrel,
                              const float &upperThresholdBarrel, const float &lowerThresholdEndcap) {
    auto lambda = [etaBoundary, lowerThresholdBarrel, upperThresholdBarrel,
                   lowerThresholdEndcap](const ROOT::RVec<float> &eta) {
        ROOT::RVec<int> mask =
            ( ( (abs(eta) >= lowerThresholdBarrel) && (abs(eta) < upperThresholdBarrel) ) || 
              ( ( abs(eta) >= lowerThresholdEndcap ) && ( abs(eta) < etaBoundary ) ) );
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {etaColumnName});
    return df1;
}
///
/// function  to make a flag that if exist dimuon pair from Higgs
ROOT::RDF::RNode DiMuonFromHiggs(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &dimuons_index) {
    auto HiggsCand_Flag = [](const ROOT::RVec<int> &dimuons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dimuons_index.at(0) == -1 || dimuons_index.at(1) == -1 ) {
                                    return 0;
                                 } else { // else exist a dimuon pair that may from Higgs
                                    return 1;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, HiggsCand_Flag, {dimuons_index});
    return df1;
}
/// function to pick dimuon pair from Higgs
ROOT::RDF::RNode HiggsToDiMuonPairCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &dimuons_index) {
    auto dimuon_calc_p4byPt = [](const ROOT::RVec<float> &particle_pts,
                                 const ROOT::RVec<float> &particle_etas,
                                 const ROOT::RVec<float> &particle_phis,
                                 const ROOT::RVec<float> &particle_masses,
                                 const ROOT::RVec<int> &dimuons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dimuons_index.at(0) == -1 || dimuons_index.at(1) == -1 ) {
                                    return ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float,default_float);
                                 } else {
                                    // int leading    mu1 = dimuons_index[0];
                                    // int subleading mu2 = dimuons_index[1];
                                    for (unsigned int k = 0; k < (int)dimuons_index.size(); ++k) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(dimuons_index[k]),
                                                                                particle_etas.at(dimuons_index[k]),
                                                                                particle_phis.at(dimuons_index[k]),
                                                                                particle_masses.at(dimuons_index[k])));
                                    }
                                    p4_dimuon = ROOT::Math::PtEtaPhiMVector((p4[0]+p4[1]).pt(),
                                                                            (p4[0]+p4[1]).eta(),
                                                                            (p4[0]+p4[1]).phi(),
                                                                            (p4[0]+p4[1]).mass());
                                    return p4_dimuon;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, dimuon_calc_p4byPt, {particle_pts, particle_etas, particle_phis, particle_masses, dimuons_index});
    return df1;
}
///
///
ROOT::RDF::RNode DiMuonFromZVeto(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &particle_charges,
                                 const std::string &goodmuons_index) {
    auto pair_calc_mass = [](const ROOT::RVec<float> &particle_pts,
                               const ROOT::RVec<float> &particle_etas,
                               const ROOT::RVec<float> &particle_phis,
                               const ROOT::RVec<float> &particle_masses,
                               const ROOT::RVec<int> &particle_charges,
                               const ROOT::RVec<int> &goodmuons_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    try {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodmuons_index[k]),
                                                                         particle_etas.at(goodmuons_index[k]),
                                                                         particle_phis.at(goodmuons_index[k]),
                                                                         particle_masses.at(goodmuons_index[k])));
                                    } catch (const std::out_of_range &e) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_1;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_2;
                                 p4_1 = p4;
                                 p4_2 = p4;
                                 for (unsigned int i = 0; i < p4_1.size(); ++i) {
                                     for (unsigned int j = i + 1; j < p4_2.size(); ++j) {
                                         if (p4_1[i].pt() < 0.0 || p4_2[j].pt() < 0.0) {
                                             continue; 
                                         }
                                         if ( particle_charges[goodmuons_index[i]] + particle_charges[goodmuons_index[j]] != 0 ) {
                                             continue;
                                         }
                                         float dimuon_mass = (p4_1[i] + p4_2[j]).mass();
                                         if (dimuon_mass >= 81.0 && dimuon_mass <= 101.0)
                                             return 0;
                                     }
                                 }
                                 return 1;
                             };
    auto df1 = 
        df.Define(outputname, pair_calc_mass, {particle_pts, particle_etas, particle_phis, particle_masses, particle_charges, goodmuons_index});
    return df1;
}
///
///
ROOT::RDF::RNode LeptonChargeSum(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_charge,
                                 const std::string &goodmuons_index) {
    auto calc_charge_sum = [](const ROOT::RVec<int> &muon_charges, 
                              const ROOT::RVec<int> &goodmuons_index) {
        int charge_sum = 0;
        // for (unsigned int i = 0; i < ele_charges.size(); ++i) {
        //     if (!std::isnan(ele_charges[i])) {
        //         charge_sum += ele_charges[i];
        //     }
        // }
        for (unsigned int i = 0; i < (int)goodmuons_index.size(); ++i) {
            if (!std::isnan(muon_charges[goodmuons_index[i]])) {
                charge_sum += muon_charges[goodmuons_index[i]];
            }
        }
        // for (unsigned int i = 0; i < tau_charges.size(); ++i) {
        //     if (!std::isnan(tau_charges[i])) {
        //         charge_sum += tau_charges[i];
        //     }
        // }
        // std::cout << "charge_sum :" << charge_sum << std::endl;
        if (charge_sum == 1 || charge_sum == -1) {
            return 1;
        } else if (charge_sum == 0) {
            return 2;
        } else {
            return 0;
        }
    };
    auto df1 =
        df.Define(outputname, calc_charge_sum, {muon_charge, goodmuons_index});
    return df1;
}
///
///
ROOT::RDF::RNode LeptonChargeSumEleMu(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_charge,
                                 const std::string &ele_charge,
                                 const std::string &goodmuons_index,
                                 const std::string &base_electrons_index) {
    auto calc_charge_sum = [](const ROOT::RVec<int> &muon_charges,
                              const ROOT::RVec<int> &ele_charges,
                              const ROOT::RVec<int> &goodmuons_index,
                              const ROOT::RVec<int> &base_electrons_index) {
        int charge_sum = 0;
        for (unsigned int i = 0; i < (int)base_electrons_index.size(); ++i) {
            if (!std::isnan(ele_charges[base_electrons_index[i]])) {
                charge_sum += ele_charges[base_electrons_index[i]];
            }
        }
        for (unsigned int i = 0; i < (int)goodmuons_index.size(); ++i) {
            if (!std::isnan(muon_charges[goodmuons_index[i]])) {
                charge_sum += muon_charges[goodmuons_index[i]];
            }
        }
        // for (unsigned int i = 0; i < tau_charges.size(); ++i) {
        //     if (!std::isnan(tau_charges[i])) {
        //         charge_sum += tau_charges[i];
        //     }
        // }
        // std::cout << "charge_sum :" << charge_sum << std::endl;
        if (charge_sum == 1 || charge_sum == -1) {
            return 1;
        } else if (charge_sum == 0) {
            return 2;
        } else {
            return 0;
        }
    };
    auto df1 =
        df.Define(outputname, calc_charge_sum, {muon_charge, ele_charge, goodmuons_index, base_electrons_index});
    return df1;
}
///
/// m2m and 4m channel to veto electrons
///
ROOT::RDF::RNode Ele_Veto(ROOT::RDF::RNode df, 
                    const std::string& output_name, 
                    const std::string& base_ele_mask) {
    auto veto_electrons = [](const ROOT::RVec<int>& ele_mask) {
        if ((int)ROOT::VecOps::Nonzero(ele_mask).size() >= 1) {
            return 0;
        } else {
            return 1;
        }
    };
    auto df1 = df.Define(output_name, veto_electrons, {base_ele_mask});
    return df1;
}
///
///
/// function to pick dimuon pair from Higgs
ROOT::RDF::RNode HiggsCandDiMuonPairCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &particle_charges,
                                 const std::string &goodmuons_index) {
    auto pair_calc_p4byPt = [](const ROOT::RVec<float> &particle_pts,
                               const ROOT::RVec<float> &particle_etas,
                               const ROOT::RVec<float> &particle_phis,
                               const ROOT::RVec<float> &particle_masses,
                               const ROOT::RVec<int> &particle_charges,
                               const ROOT::RVec<int> &goodmuons_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    try {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodmuons_index[k]),  ///goodmuons_index[k] points to the good muon index k
                                                                         particle_etas.at(goodmuons_index[k]),          // k = 0, points to goodmuon_index[0]
                                                                         particle_phis.at(goodmuons_index[k]),          // k ,points to goodmuon_index[k]
                                                                         particle_masses.at(goodmuons_index[k])));      // index what I want is goodmuon_index[k] k,i or j
                                    } catch (const std::out_of_range &e) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_1;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_2;
                                 p4_1 = p4;
                                 p4_2 = p4;
                                 float ptsum = -1;
                                 int index1 = -1,index2 = -1;
                                 for (unsigned int i = 0; i < p4_1.size(); ++i) {
                                     for (unsigned int j = i + 1; j < p4_2.size(); ++j) {
                                         if (p4_1[i].pt() < 0.0 || p4_2[j].pt() < 0.0)
                                             continue; 
                                         /// need opposite sign dimuons
                                         if ( particle_charges[goodmuons_index[i]] + particle_charges[goodmuons_index[j]] != 0 ) {
                                             continue;
                                         }
                                         /// Add dimuon mass window
                                         if ( (p4_1[i] + p4_2[j]).mass() < 110 || (p4_1[i] + p4_2[j]).mass() > 150 ) {
                                             continue;
                                         }
                                         if ( p4_1[i].pt() + p4_2[j].pt() > ptsum) {
                                             ptsum = p4_1[i].pt() + p4_2[j].pt();
                                             if ( p4_1[i].pt() > p4_1[j].pt() ) {
                                                index1 = goodmuons_index[i];
                                                index2 = goodmuons_index[j];
                                             } else {
                                                index1 = goodmuons_index[j];
                                                index2 = goodmuons_index[i]; /// need to return the index1 and index2 as goodmuons pair collection.
                                             }
                                         }
                                     }
                                 }
                                 ROOT::RVec<int> DiMuonPair = {index1, index2};
                                 return DiMuonPair;
                                 ///p4_dimuon = p4_dileptonsystem[0];
                                 ///return p4_dimuon; /// return dimuon_pair_p4 order by pt
                             };
    auto df1 = 
        df.Define(outputname, pair_calc_p4byPt, {particle_pts, particle_etas, particle_phis, particle_masses, particle_charges, goodmuons_index});
    return df1;
}
///
/// need Zee mass cut in [70,110]
ROOT::RDF::RNode ZCandDiElectronPairCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &particle_charges,
                                 const std::string &base_electrons_index) {
    auto pair_calc_p4byPt = [](const ROOT::RVec<float> &particle_pts,
                               const ROOT::RVec<float> &particle_etas,
                               const ROOT::RVec<float> &particle_phis,
                               const ROOT::RVec<float> &particle_masses,
                               const ROOT::RVec<int> &particle_charges,
                               const ROOT::RVec<int> &base_electrons_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 for (unsigned int k = 0; k < (int)base_electrons_index.size(); ++k) {
                                    try {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(base_electrons_index[k]),
                                                                         particle_etas.at(base_electrons_index[k]),         
                                                                         particle_phis.at(base_electrons_index[k]),          
                                                                         particle_masses.at(base_electrons_index[k])));      
                                    } catch (const std::out_of_range &e) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_1;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_2;
                                 p4_1 = p4;
                                 p4_2 = p4;
                                 float ptsum = -1;
                                 int index1 = -1,index2 = -1;
                                 for (unsigned int i = 0; i < p4_1.size(); ++i) {
                                     for (unsigned int j = i + 1; j < p4_2.size(); ++j) {
                                         if (p4_1[i].pt() < 0.0 || p4_2[j].pt() < 0.0)
                                             continue; 
                                         /// need opposite sign dimuons
                                         if ( particle_charges[base_electrons_index[i]] + particle_charges[base_electrons_index[j]] != 0 ) {
                                             continue;
                                         }
                                         /// Add dimuon mass window
                                         if ( (p4_1[i] + p4_2[j]).mass() < 70 || (p4_1[i] + p4_2[j]).mass() > 110 ) {
                                             continue;
                                         }
                                         if ( p4_1[i].pt() + p4_2[j].pt() > ptsum) {
                                             ptsum = p4_1[i].pt() + p4_2[j].pt();
                                             if ( p4_1[i].pt() > p4_1[j].pt() ) {
                                                index1 = base_electrons_index[i];
                                                index2 = base_electrons_index[j];
                                             } else {
                                                index1 = base_electrons_index[j];
                                                index2 = base_electrons_index[i]; /// need to return the index1 and index2 as goodmuons pair collection.
                                             }
                                         }
                                     }
                                 }
                                 ROOT::RVec<int> DiElectronPair = {index1, index2};
                                 return DiElectronPair;
                             };
    auto df1 = 
        df.Define(outputname, pair_calc_p4byPt, {particle_pts, particle_etas, particle_phis, particle_masses, particle_charges, base_electrons_index});
    return df1;
}
/// function to make a flag that if exist dielectron pair from Z
ROOT::RDF::RNode DiEleFromZ(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &dielectrons_index) {
    auto ZCand_Flag = [](const ROOT::RVec<int> &dielectrons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dielectrons_index.at(0) == -1 || dielectrons_index.at(1) == -1 ) {
                                    return 0;
                                 } else { // else exist a diele pair that may from Z
                                    return 1;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, ZCand_Flag, {dielectrons_index});
    return df1;
}
/// function to pick dimuon pair from Higgs
ROOT::RDF::RNode ZToDiElectronPairCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &dielectrons_index) {
    auto dielectron_calc_p4byPt = [](const ROOT::RVec<float> &particle_pts,
                                 const ROOT::RVec<float> &particle_etas,
                                 const ROOT::RVec<float> &particle_phis,
                                 const ROOT::RVec<float> &particle_masses,
                                 const ROOT::RVec<int> &dielectrons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_diele;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dielectrons_index.at(0) == -1 || dielectrons_index.at(1) == -1 ) {
                                    return ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float,default_float);
                                 } else {
                                    // int leading    ele1 = dielectrons_index[0];
                                    // int subleading ele2 = dielectrons_index[1];
                                    for (unsigned int k = 0; k < (int)dielectrons_index.size(); ++k) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(dielectrons_index[k]),
                                                                                particle_etas.at(dielectrons_index[k]),
                                                                                particle_phis.at(dielectrons_index[k]),
                                                                                particle_masses.at(dielectrons_index[k])));
                                    }
                                    p4_diele = ROOT::Math::PtEtaPhiMVector((p4[0]+p4[1]).pt(),
                                                                            (p4[0]+p4[1]).eta(),
                                                                            (p4[0]+p4[1]).phi(),
                                                                            (p4[0]+p4[1]).mass());
                                    return p4_diele;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, dielectron_calc_p4byPt, {particle_pts, particle_etas, particle_phis, particle_masses, dielectrons_index});
    return df1;
}
///
/// \\\NOTICE
/// 4m events has only 4 muons
///
/// In 4m events, if it is possible to form two distinct m+m- pair 
/// each with a mass between 81 and 101 GeV, the event is discarded.
/// 
/// In 4m events, one muon pair must have mass between 110 and 150 GeV, 
/// and the other muon pair must have mass between 81 and 101 GeV.
///
/// In 4m events, if both combinations have a muon pair in the Z-mass window 
/// and a muon pair in the signal-mass window, 
/// the combination in which the mass of the Z candidate is closer to 91 GeV is chosen.
///
/// This function is need to make a mask (return the muon index)
/// make two flag. 
/// notice that 1,2,3,4 muons can make 6 types of pair. 
/// notice if SFOS, only 2 types of pair
ROOT::RDF::RNode HiggsAndZFourMuonsCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &particle_charges,
                                 const std::string &goodmuons_index) {
    auto pair_calc_p4 = [](const ROOT::RVec<float> &particle_pts,
                               const ROOT::RVec<float> &particle_etas,
                               const ROOT::RVec<float> &particle_phis,
                               const ROOT::RVec<float> &particle_masses,
                               const ROOT::RVec<int> &particle_charges,
                               const ROOT::RVec<int> &goodmuons_index) {
                                std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    try {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodmuons_index[k]),  ///goodmuons_index[k] points to the good muon index k
                                                                         particle_etas.at(goodmuons_index[k]),          // k = 0, points to goodmuon_index[0]
                                                                         particle_phis.at(goodmuons_index[k]),          // k ,points to goodmuon_index[k]
                                                                         particle_masses.at(goodmuons_index[k])));      // index what I want is goodmuon_index[k] k,i or j
                                    } catch (const std::out_of_range &e) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                }
                                int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
                                float min_mass_diff = 999999.0;
                                for (unsigned int i = 0; i < p4.size(); ++i) {
                                    for (unsigned int j = i + 1; j < p4.size(); ++j) {
                                        int pair1_index1 = i;
                                        int pair1_index2 = j;
                                        int pair2_index1 = -1, pair2_index2 = -1;
                                        for (unsigned int k = 0; k < p4.size(); ++k) {
                                            if (k != i && k != j) {
                                                if (pair2_index1 == -1) {
                                                    pair2_index1 = k;
                                                } else {
                                                    pair2_index2 = k;
                                                }
                                            }
                                        }
                                        if (p4[pair1_index1].pt() < 0.0 || p4[pair1_index2].pt() < 0.0 
                                        || p4[pair2_index1].pt() < 0.0 || p4[pair2_index2].pt() < 0.0) {
                                            continue;
                                        }
                                        /// need SFOS
                                        if ( particle_charges[goodmuons_index[pair1_index1]] + particle_charges[goodmuons_index[pair1_index2]] != 0 ) {
                                            continue;
                                        }
                                        if ( particle_charges[goodmuons_index[pair2_index1]] + particle_charges[goodmuons_index[pair2_index2]] != 0 ) {
                                            continue;
                                        }
                                        /// two muon mass and the other two muon
                                        float dimuon_mass_1 = (p4[pair1_index1] + p4[pair1_index2]).mass();
                                        float dimuon_mass_2 = (p4[pair2_index1] + p4[pair2_index2]).mass();
                                        ///
                                        /// Another, if two pairs both in Z window, return {-1,-1,-1,-1} at once
                                        if ( dimuon_mass_1 > 81 && dimuon_mass_1 < 101 ) {
                                            if (dimuon_mass_2 > 81 && dimuon_mass_2 < 101) {
                                                return ROOT::RVec<int> {-1,-1,-1,-1};
                                            }
                                        }
                                        ///
                                        /// this if can select one pair in Higgs and another in Z
                                        if ( dimuon_mass_1 > 110 && dimuon_mass_1 < 150 ) {
                                            if (dimuon_mass_2 > 81 && dimuon_mass_2 < 101) {
                                                /// need Z cand pair close to 91
                                                /// the first time will always be true in below if
                                                /// this if can select the ZCand which dimuon mass closer to 91 GeV
                                                if ( fabs(dimuon_mass_2 - 91) < fabs( min_mass_diff -91 ) ) {
                                                    min_mass_diff = dimuon_mass_2;
                                                    /// pt ordering
                                                    if ( p4[pair1_index1].pt() > p4[pair1_index2].pt() ) {
                                                        index1 = goodmuons_index[pair1_index1];
                                                        index2 = goodmuons_index[pair1_index2];
                                                    } else {
                                                        index1 = goodmuons_index[pair1_index2];
                                                        index2 = goodmuons_index[pair1_index1];
                                                    }
                                                    if ( p4[pair2_index1].pt() > p4[pair2_index2].pt() ) {
                                                        index3 = goodmuons_index[pair2_index1];
                                                        index4 = goodmuons_index[pair2_index2];
                                                    } else {
                                                        index3 = goodmuons_index[pair2_index2];
                                                        index4 = goodmuons_index[pair2_index1];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                /// if no two pairs both in Z window among all the loop, return the four goodmuons index
                                ROOT::RVec<int> fourmuons_idx = {index1, index2, index3, index4};
                                return fourmuons_idx;
                            };
    auto df1 = 
        df.Define(outputname, pair_calc_p4, {particle_pts, particle_etas, particle_phis, particle_masses, particle_charges, goodmuons_index});
    return df1;
}
///
///
ROOT::RDF::RNode QuadMuonFromZZVeto(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &quadmuons_index) {
    auto HiggsZCand_Flag = [](const ROOT::RVec<int> &quadmuons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( quadmuons_index.at(0) == -1 || quadmuons_index.at(1) == -1 || quadmuons_index.at(2) == -1 || quadmuons_index.at(3) == -1 ) {
                                    return 0;
                                 } else { // else exist four muons that may from Higgs and Z
                                    return 1;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, HiggsZCand_Flag, {quadmuons_index});
    return df1;
}
///
/// pick the third and fourth muon as Z Cand
ROOT::RDF::RNode ZToSecondMuonPairCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &fourmuons_index) {
    auto dimuonZ_calc_p4 = [](const ROOT::RVec<float> &particle_pts,
                                 const ROOT::RVec<float> &particle_etas,
                                 const ROOT::RVec<float> &particle_phis,
                                 const ROOT::RVec<float> &particle_masses,
                                 const ROOT::RVec<int> &fourmuons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuZ;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( fourmuons_index.at(2) == -1 || fourmuons_index.at(3) == -1 ) {
                                    return ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float,default_float);
                                 } else {
                                    // int leading    Zmu1 = fourmuons_index[2];
                                    // int subleading Zmu2 = fourmuons_index[3];
                                    for (unsigned int k = 2; k < (int)fourmuons_index.size(); ++k) {
                                        /// push back [2] and [3] 
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(fourmuons_index[k]),
                                                                                particle_etas.at(fourmuons_index[k]),
                                                                                particle_phis.at(fourmuons_index[k]),
                                                                                particle_masses.at(fourmuons_index[k])));
                                    }
                                    p4_dimuZ = ROOT::Math::PtEtaPhiMVector((p4[0]+p4[1]).pt(),
                                                                            (p4[0]+p4[1]).eta(),
                                                                            (p4[0]+p4[1]).phi(),
                                                                            (p4[0]+p4[1]).mass());
                                    return p4_dimuZ;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, dimuonZ_calc_p4, {particle_pts, particle_etas, particle_phis, particle_masses, fourmuons_index});
    return df1;
}
///
/// extra muon in the m2m channel
ROOT::RDF::RNode ExtraMuonIndexFromW(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &goodmuons_index,
                                 const std::string &dimuons_index) {
    auto extra_muonW_calc_p4 = [](const ROOT::RVec<float> &particle_pts,
                                const ROOT::RVec<float> &particle_etas,
                                const ROOT::RVec<float> &particle_phis,
                                const ROOT::RVec<float> &particle_masses,
                                const ROOT::RVec<int> &goodmuons_index,
                                const ROOT::RVec<int> &dimuons_index) {
                                ROOT::Math::PtEtaPhiMVector p4_dimuZ;
                                std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                /// dimuons_index [0] and [1] stands the Higgs Mu1 and Mu2
                                /// goodmuons_index [0] [1] [2] stands 3 muons
                                /// do a loop to find the extra muon index 
                                for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    if ( goodmuons_index[k] != dimuons_index[0] && goodmuons_index[k] != dimuons_index[1] ) {
                                        if ( dimuons_index[0] != -1 && dimuons_index[1] != -1 ) {
                                            // vec(1,6); means "1" stands length, "6" stands inital value
                                            ROOT::RVec<int> extra_muon_index(1, goodmuons_index[k]);
                                            return extra_muon_index;
                                        }
                                    }
                                }
                                ROOT::RVec<int> NO_Index(1, -1);
                                return NO_Index;
                            };
    auto df1 = 
        df.Define(outputname, extra_muonW_calc_p4, {particle_pts, particle_etas, particle_phis, particle_masses, goodmuons_index, dimuons_index});
    return df1;
}
///
/// extra muon p4 in the m2m channel
ROOT::RDF::RNode ExtraMuonFromW(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &extra_muon_index) {
    auto extra_muonW_calc_p4 = [](const ROOT::RVec<float> &particle_pts,
                                const ROOT::RVec<float> &particle_etas,
                                const ROOT::RVec<float> &particle_phis,
                                const ROOT::RVec<float> &particle_masses,
                                const ROOT::RVec<int> &extra_muon_index) {
                                std::vector<ROOT::Math::PtEtaPhiMVector> p4_lepFromW;
                                if ( extra_muon_index[0] == -1 ) {
                                    p4_lepFromW.push_back( ROOT::Math::PtEtaPhiMVector(default_float,default_float,default_float,default_float) );
                                } else {
                                    p4_lepFromW.push_back( ROOT::Math::PtEtaPhiMVector(particle_pts.at(extra_muon_index[0]),
                                                                                    particle_etas.at(extra_muon_index[0]),
                                                                                    particle_phis.at(extra_muon_index[0]),
                                                                                    particle_masses.at(extra_muon_index[0])) );
                                }
                                return p4_lepFromW[0];
                            };
    auto df1 = 
        df.Define(outputname, extra_muonW_calc_p4, {particle_pts, particle_etas, particle_phis, particle_masses, extra_muon_index});
    return df1;
}
/// function that calc the p4 of the SS or OS mu compared to extra lepton from W
ROOT::RDF::RNode muSSorOSwithLeptonW_p4(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_pts,
                                 const std::string &muon_etas,
                                 const std::string &muon_phis,
                                 const std::string &muon_masses,
                                 const std::string &muon_charges,
                                 const std::string &lep_charges,
                                 const std::string &dimuons_index,
                                 const std::string &lep_index,
                                 const int SameSign) {
    auto calc_p4 = [SameSign](const ROOT::RVec<float> &muon_pts,
                                 const ROOT::RVec<float> &muon_etas,
                                 const ROOT::RVec<float> &muon_phis,
                                 const ROOT::RVec<float> &muon_masses,
                                 const ROOT::RVec<int> &muon_charges,
                                 const ROOT::RVec<int> &lep_charges,
                                 const ROOT::RVec<int> &dimuons_index,
                                 const ROOT::RVec<int> &lep_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dimuons_index.at(0) == -1 || dimuons_index.at(1) == -1 || lep_index.at(0) == -1 ) {
                                    return ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float,default_float);
                                 } else {
                                    if ( SameSign == 1 ) {
                                        for (unsigned int k = 0; k < (int)dimuons_index.size(); ++k) {
                                            /// same sign
                                            if ( lep_charges[lep_index[0]] ==  muon_charges[dimuons_index[k]] ) {
                                                p4.push_back(ROOT::Math::PtEtaPhiMVector(muon_pts.at(dimuons_index[k]),
                                                                                        muon_etas.at(dimuons_index[k]),
                                                                                        muon_phis.at(dimuons_index[k]),
                                                                                        muon_masses.at(dimuons_index[k])));
                                            }
                                        }
                                    } else {
                                        for (unsigned int k = 0; k < (int)dimuons_index.size(); ++k) {
                                            /// opposite sign
                                            if ( lep_charges[lep_index[0]] + muon_charges[dimuons_index[k]] == 0) {
                                                p4.push_back(ROOT::Math::PtEtaPhiMVector(muon_pts.at(dimuons_index[k]),
                                                                                        muon_etas.at(dimuons_index[k]),
                                                                                        muon_phis.at(dimuons_index[k]),
                                                                                        muon_masses.at(dimuons_index[k])));
                                            }
                                        }
                                    }
                                    return p4[0];
                                 }
                             };
    auto df1 = 
        df.Define(outputname, calc_p4, {muon_pts, muon_etas, muon_phis, muon_masses, muon_charges, lep_charges, dimuons_index, lep_index});
    return df1;
}
///
///// function to calc the MHT 
ROOT::RDF::RNode MHT_Calculation(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &goodjets_index) {
    auto MHT_calc_p4 = [](const ROOT::RVec<float> &particle_pts,
                                 const ROOT::RVec<float> &particle_etas,
                                 const ROOT::RVec<float> &particle_phis,
                                 const ROOT::RVec<float> &particle_masses,
                                 const ROOT::RVec<int> &goodjets_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_MHT(0, 0, 0, 0);
                                 ///std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 ///
                                //  for (unsigned int i = 0; i < (int)goodjets_index.size(); ++i ) {
                                //     try {
                                //         p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodjets_index[i]), 
                                //                                          particle_etas.at(goodjets_index[i]),
                                //                                          particle_phis.at(goodjets_index[i]),
                                //                                          particle_masses.at(goodjets_index[i])));
                                //     } catch (const std::out_of_range &e) {
                                //         p4.push_back(ROOT::Math::PtEtaPhiMVector(0,0,0,0));
                                //     }
                                //  }
                                 for (unsigned int j = 0; j < (int)goodjets_index.size(); ++j ) {
                                    if ( particle_pts.at(goodjets_index[j]) > 30 && fabs( particle_etas.at(goodjets_index[j]) ) < 4.7 ) {
                                        p4_MHT += ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodjets_index[j]), 
                                                                         particle_etas.at(goodjets_index[j]),
                                                                         particle_phis.at(goodjets_index[j]),
                                                                         particle_masses.at(goodjets_index[j]));
                                    }
                                 }
                                 ///
                                 p4_MHT = -p4_MHT;
                                 ///p4_MHT = -ROOT::Math::PtEtaPhiEVector(p4_MHT.Pt(), 0, p4_MHT.Phi(), p4_MHT.Pt());
                                 return (ROOT::Math::PtEtaPhiMVector)p4_MHT;
                             };
    auto df1 = 
        df.Define(outputname, MHT_calc_p4, {particle_pts, particle_etas, particle_phis, particle_masses, goodjets_index});
    return df1;
}
/// UF calculation about MHT?
ROOT::RDF::RNode MHT_CalculationALL(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_pts,
                                 const std::string &muon_etas,
                                 const std::string &muon_phis,
                                 const std::string &muon_masses,
                                 const std::string &muon_index,
                                 const std::string &ele_pts,
                                 const std::string &ele_etas,
                                 const std::string &ele_phis,
                                 const std::string &ele_masses,
                                 const std::string &ele_index,
                                 const std::string &jet_pts,
                                 const std::string &jet_etas,
                                 const std::string &jet_phis,
                                 const std::string &jet_masses,
                                 const std::string &goodjets_index) {
    auto MHT_calc_p4 = [](const ROOT::RVec<float> &muon_pts,
                                 const ROOT::RVec<float> &muon_etas,
                                 const ROOT::RVec<float> &muon_phis,
                                 const ROOT::RVec<float> &muon_masses,
                                 const ROOT::RVec<int> &muon_index,
                                 const ROOT::RVec<float> &ele_pts,
                                 const ROOT::RVec<float> &ele_etas,
                                 const ROOT::RVec<float> &ele_phis,
                                 const ROOT::RVec<float> &ele_masses,
                                 const ROOT::RVec<int> &ele_index,
                                 const ROOT::RVec<float> &jet_pts,
                                 const ROOT::RVec<float> &jet_etas,
                                 const ROOT::RVec<float> &jet_phis,
                                 const ROOT::RVec<float> &jet_masses,
                                 const ROOT::RVec<int> &goodjets_index) {
                                ROOT::Math::PtEtaPhiMVector p4_MHTALL(0, 0, 0, 0);
                                ///std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                /// muon
                                for (unsigned int j = 0; j < (int)muon_index.size(); ++j ) {
                                        p4_MHTALL += ROOT::Math::PtEtaPhiMVector(muon_pts.at(muon_index[j]), 
                                                                        muon_etas.at(muon_index[j]),
                                                                        muon_phis.at(muon_index[j]),
                                                                        muon_masses.at(muon_index[j]));
                                }
                                /// ele
                                for (unsigned int j = 0; j < (int)ele_index.size(); ++j ) {
                                        p4_MHTALL += ROOT::Math::PtEtaPhiMVector(ele_pts.at(ele_index[j]), 
                                                                        ele_etas.at(ele_index[j]),
                                                                        ele_phis.at(ele_index[j]),
                                                                        ele_masses.at(ele_index[j]));
                                }
                                /// jet
                                for (unsigned int j = 0; j < (int)goodjets_index.size(); ++j ) {
                                    if ( jet_pts.at(goodjets_index[j]) > 30 && fabs( jet_etas.at(goodjets_index[j]) ) < 4.7 ) {
                                        p4_MHTALL += ROOT::Math::PtEtaPhiMVector(jet_pts.at(goodjets_index[j]), 
                                                                        jet_etas.at(goodjets_index[j]),
                                                                        jet_phis.at(goodjets_index[j]),
                                                                        jet_masses.at(goodjets_index[j]));
                                    }
                                }
                                ///
                                p4_MHTALL = -p4_MHTALL;
                                ///p4_MHT = -ROOT::Math::PtEtaPhiEVector(p4_MHT.Pt(), 0, p4_MHT.Phi(), p4_MHT.Pt());
                                return (ROOT::Math::PtEtaPhiMVector)p4_MHTALL;
                             };
    auto df1 = 
        df.Define(outputname, MHT_calc_p4, {muon_pts, muon_etas, muon_phis, muon_masses, muon_index, ele_pts, ele_etas, ele_phis, ele_masses, ele_index, jet_pts, jet_etas, jet_phis, jet_masses, goodjets_index});
    return df1;
}
///
/// function to make all pass flagZmassVeto in e2m that be the same as m2m
ROOT::RDF::RNode PassFlag(ROOT::RDF::RNode df, const std::string &outputname) {
    auto Pass_Flag = []() {
                            return 1;
                        };
    auto df1 = 
        df.Define(outputname, Pass_Flag, {});
    return df1;
}
///
ROOT::RDF::RNode PassDiEleIn4m(ROOT::RDF::RNode df, const std::string &outputname) {
    auto Pass = []() {
                            return 999.0f;
                        };
    auto df1 = 
        df.Define(outputname, Pass, {});
    return df1;
}
///
ROOT::RDF::RNode Calc_CosThetaStar(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &lepton_p4, const std::string &mu_p4) {
    auto calculate_costhstar = [](ROOT::Math::PtEtaPhiMVector &lep_p4,
                           ROOT::Math::PtEtaPhiMVector &mu_p4) {
        TLorentzVector lepton_p4_TL;
        TLorentzVector mu_p4_TL;

        lepton_p4_TL.SetPtEtaPhiM(lep_p4.Pt(), lep_p4.Eta(), lep_p4.Phi(), lep_p4.M());
        mu_p4_TL.SetPtEtaPhiM(mu_p4.Pt(), mu_p4.Eta(), mu_p4.Phi(), mu_p4.M());

        TLorentzVector lepmu = lepton_p4_TL + mu_p4_TL;


        TVector3 lepmu_v = lepmu.Vect();
        TVector3 WHboost = -(lepmu.BoostVector());
        lepton_p4_TL.Boost(WHboost);
        mu_p4_TL.Boost(WHboost);
        TVector3 lep_v = lepton_p4_TL.Vect();
        TVector3 muos_v = mu_p4_TL.Vect();

        float cosh_angle = cos(lep_v.Angle(lepmu_v));

        if ( !std::isnan(cosh_angle) && !std::isinf(cosh_angle) ) {
            return cosh_angle;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calculate_costhstar, {lepton_p4, mu_p4});
}
///
///
ROOT::RDF::RNode Calc_CosThetaStar_ZH(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &Z_p4, const std::string &H_p4) {
    auto calculate_costhstar = [](ROOT::Math::PtEtaPhiMVector &Z_p4,
                           ROOT::Math::PtEtaPhiMVector &H_p4) {
        TLorentzVector Z_p4_TL;
        TLorentzVector H_p4_TL;

        Z_p4_TL.SetPtEtaPhiM(Z_p4.Pt(), Z_p4.Eta(), Z_p4.Phi(), Z_p4.M());
        H_p4_TL.SetPtEtaPhiM(H_p4.Pt(), H_p4.Eta(), H_p4.Phi(), H_p4.M());
        TLorentzVector TL = Z_p4_TL +H_p4_TL;

        TVector3 ZH_v = TL.Vect();
        TVector3 ZHboost = -(TL.BoostVector());
        Z_p4_TL.Boost(ZHboost);
        H_p4_TL.Boost(ZHboost);
        TVector3 Z_v = Z_p4_TL.Vect();
        TVector3 H_v = H_p4_TL.Vect();

        float cosh_angle = cos(Z_v.Angle(ZH_v));

        if ( !std::isnan(cosh_angle) && !std::isinf(cosh_angle) ) {
            return cosh_angle;
        } else {
            return -10.0f;
        }
    };
    return df.Define(outputname, calculate_costhstar, {Z_p4, H_p4});
}
///
ROOT::RDF::RNode RedirectZlepID(ROOT::RDF::RNode df, const int ifMu, const std::string &outputname) {
    auto LepID = [ifMu]() {
        if (ifMu == 1) {
            return 13;
        } else if (ifMu == 0) {
            return 11;
        } else {
            return 0; 
        }
    };
    auto df1 = 
        df.Define(outputname, LepID, {});
    return df1;
}
///
//// met_p4 PtEtaPhiM
//// met_p4 only Pt and Phi
ROOT::RDF::RNode MetCut(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &met_p4, const float pt_threshold) {
    auto cut_met = [pt_threshold](ROOT::Math::PtEtaPhiMVector &met) {
                                if ((float)met.Pt() >= pt_threshold) {
                                    return 1;
                                } else {
                                    return 0;
                                }
                            };
    auto df1 = 
        df.Define(outputname, cut_met, {met_p4});
    return df1;
}
ROOT::RDF::RNode MaxMetCut(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &met_p4, const float pt_threshold) {
    auto cut_met = [pt_threshold](ROOT::Math::PtEtaPhiMVector &met) {
                                if ((float)met.Pt() <= pt_threshold) {
                                    return 1;
                                } else {
                                    return 0;
                                }
                            };
    auto df1 = 
        df.Define(outputname, cut_met, {met_p4});
    return df1;
}

////
/// function to pick dimuon Gen pair from Higgs
ROOT::RDF::RNode HiggsCandDiMuonGenPairCollection(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &Muon_indexToGen,
                                 const std::string &dimuon_index) {
    auto RecoToGen = [](const ROOT::RVec<int> &Muon_indexToGen,
                               const ROOT::RVec<int> &dimuon_index) {
                                 int index1 = -1,index2 = -1;
                                 index1 = Muon_indexToGen.at(dimuon_index[0]);
                                 index2 = Muon_indexToGen.at(dimuon_index[1]);
                                 ROOT::RVec<int> DiMuonGenPair = {index1, index2};
                                 return DiMuonGenPair; // Two gen muon index 
                             };
    auto df1 = 
        df.Define(outputname, RecoToGen, {Muon_indexToGen, dimuon_index});
    return df1;
}
///
/// function to find the Vector Boson decay
ROOT::RDF::RNode BosonDecayMode(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &GenPart_pdgId,
                                 const std::string &GenPart_motherid,
                                 const std::string &GenPart_statusFlags) {
    auto DecayMode = [](const ROOT::RVec<int> &GenPart_pdgId,
                        const ROOT::RVec<short> &GenPart_motherid,
                        const ROOT::RVec<unsigned short> &GenPart_statusFlags) {
                                 const ROOT::RVec<short> ids_lepton = {11, 12, 13, 14, 15, 16};
                                 const ROOT::RVec<short> ids_hadron = {1, 2, 3, 4, 5, 6, 21};
                                 /// #isLastCopy = ( ((part.statusFlags) >> 13) & 1 )
                                 const ROOT::RVec<short> ids_lepton_Z = {11, 13, 15};
                                 const ROOT::RVec<short> ids_hadron_Z = {1, 2, 3, 4, 5, 6, 21};
                                 const ROOT::RVec<short> ids_invis = {12, 14, 16};
                                 for (unsigned int i = 0; i < (int)GenPart_pdgId.size(); ++i ) {
                                    // check if the gen particle is W
                                    if ( ( abs(GenPart_pdgId.at(i)) == 24 ) && ( (GenPart_statusFlags.at(i) >> 13) & 1 == 1 ) ) {
                                        // check gen particle's daughter
                                        for (unsigned int j = 0; j < (int)GenPart_pdgId.size(); ++j ) {
                                            // find the gen part daughter
                                            if ( GenPart_motherid.at(j) == i ) {
                                                // return abs(GenPart_pdgId.at(j));
                                                if (std::find(ids_lepton.begin(), ids_lepton.end(), abs(GenPart_pdgId.at(j))) != ids_lepton.end()) {
                                                    return 0; // 0 stands for W leptonic decay
                                                } else if (std::find(ids_hadron.begin(), ids_hadron.end(), abs(GenPart_pdgId.at(j))) != ids_hadron.end()) {
                                                    return 1; // 1 stands for W hadronic decay
                                                }
                                            }
                                        }
                                        return -1; // -1 stands no W or Z
                                    } else if ( ( abs(GenPart_pdgId.at(i)) == 23 ) && ( (GenPart_statusFlags.at(i) >> 13) & 1 == 1 ) ) { // gen particle is Z
                                        // check gen particle's daughter
                                        for (unsigned int k = 0; k < (int)GenPart_pdgId.size(); ++k ) {
                                            // find the gen part daughter
                                            if ( GenPart_motherid.at(k) == i ) {
                                                // return abs(GenPart_pdgId.at(k));
                                                if (std::find(ids_lepton_Z.begin(), ids_lepton_Z.end(), abs(GenPart_pdgId.at(k))) != ids_lepton_Z.end()) {
                                                    return 0; // 0 stands for Z leptonic decay
                                                } else if (std::find(ids_hadron_Z.begin(), ids_hadron_Z.end(), abs(GenPart_pdgId.at(k))) != ids_hadron_Z.end()) {
                                                    return 1; // 1 stands for Z hadronic decay
                                                } else if (std::find(ids_invis.begin(), ids_invis.end(), abs(GenPart_pdgId.at(k))) != ids_invis.end()) {
                                                    return 2; // 2 stands for Z invisible decay
                                                }
                                            }
                                        }
                                        return -1; // -1 stands no W or Z
                                    } 
                                 }
                                 return -1; // -1 stands no W or Z
                             };
    auto df1 = 
        df.Define(outputname, DecayMode, {GenPart_pdgId, GenPart_motherid, GenPart_statusFlags});
    return df1;
}
ROOT::RDF::RNode BosonDecayMode_run2(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &GenPart_pdgId,
                                 const std::string &GenPart_motherid,
                                 const std::string &GenPart_statusFlags) {
    auto DecayMode = [](const ROOT::RVec<int> &GenPart_pdgId,
                        const ROOT::RVec<int> &GenPart_motherid,
                        const ROOT::RVec<int> &GenPart_statusFlags) {
                                 const ROOT::RVec<int> ids_lepton = {11, 12, 13, 14, 15, 16};
                                 const ROOT::RVec<int> ids_hadron = {1, 2, 3, 4, 5, 6, 21};
                                 /// #isLastCopy = ( ((part.statusFlags) >> 13) & 1 )
                                 const ROOT::RVec<int> ids_lepton_Z = {11, 13, 15};
                                 const ROOT::RVec<int> ids_hadron_Z = {1, 2, 3, 4, 5, 6, 21};
                                 const ROOT::RVec<int> ids_invis = {12, 14, 16};
                                 for (unsigned int i = 0; i < (int)GenPart_pdgId.size(); ++i ) {
                                    // check if the gen particle is W
                                    if ( ( abs(GenPart_pdgId.at(i)) == 24 ) && ( (GenPart_statusFlags.at(i) >> 13) & 1 == 1 ) ) {
                                        // check gen particle's daughter
                                        for (unsigned int j = 0; j < (int)GenPart_pdgId.size(); ++j ) {
                                            // find the gen part daughter
                                            if ( GenPart_motherid.at(j) == i ) {
                                                // return abs(GenPart_pdgId.at(j));
                                                if (std::find(ids_lepton.begin(), ids_lepton.end(), abs(GenPart_pdgId.at(j))) != ids_lepton.end()) {
                                                    return 0; // 0 stands for W leptonic decay
                                                } else if (std::find(ids_hadron.begin(), ids_hadron.end(), abs(GenPart_pdgId.at(j))) != ids_hadron.end()) {
                                                    return 1; // 1 stands for W hadronic decay
                                                }
                                            }
                                        }
                                        return -1; // -1 stands no W or Z
                                    } else if ( ( abs(GenPart_pdgId.at(i)) == 23 ) && ( (GenPart_statusFlags.at(i) >> 13) & 1 == 1 ) ) { // gen particle is Z
                                        // check gen particle's daughter
                                        for (unsigned int k = 0; k < (int)GenPart_pdgId.size(); ++k ) {
                                            // find the gen part daughter
                                            if ( GenPart_motherid.at(k) == i ) {
                                                // return abs(GenPart_pdgId.at(k));
                                                if (std::find(ids_lepton_Z.begin(), ids_lepton_Z.end(), abs(GenPart_pdgId.at(k))) != ids_lepton_Z.end()) {
                                                    return 0; // 0 stands for Z leptonic decay
                                                } else if (std::find(ids_hadron_Z.begin(), ids_hadron_Z.end(), abs(GenPart_pdgId.at(k))) != ids_hadron_Z.end()) {
                                                    return 1; // 1 stands for Z hadronic decay
                                                } else if (std::find(ids_invis.begin(), ids_invis.end(), abs(GenPart_pdgId.at(k))) != ids_invis.end()) {
                                                    return 2; // 2 stands for Z invisible decay
                                                }
                                            }
                                        }
                                        return -1; // -1 stands no W or Z
                                    } 
                                 }
                                 return -1; // -1 stands no W or Z
                             };
    auto df1 = 
        df.Define(outputname, DecayMode, {GenPart_pdgId, GenPart_motherid, GenPart_statusFlags});
    return df1;
}
/// function to pick dimuon pair in DY control region
ROOT::RDF::RNode DY_DiMuonPair_CR(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &particle_charges,
                                 const std::string &goodmuons_index) {
    auto pair_calc_p4byPt = [](const ROOT::RVec<float> &particle_pts,
                               const ROOT::RVec<float> &particle_etas,
                               const ROOT::RVec<float> &particle_phis,
                               const ROOT::RVec<float> &particle_masses,
                               const ROOT::RVec<int> &particle_charges,
                               const ROOT::RVec<int> &goodmuons_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    try {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(goodmuons_index[k]),  ///goodmuons_index[k] points to the good muon index k
                                                                         particle_etas.at(goodmuons_index[k]),          // k = 0, points to goodmuon_index[0]
                                                                         particle_phis.at(goodmuons_index[k]),          // k ,points to goodmuon_index[k]
                                                                         particle_masses.at(goodmuons_index[k])));      // index what I want is goodmuon_index[k] k,i or j
                                    } catch (const std::out_of_range &e) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_1;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_2;
                                 p4_1 = p4;
                                 p4_2 = p4;
                                 float ptsum = -1;
                                 int index1 = -1,index2 = -1;
                                 for (unsigned int i = 0; i < p4_1.size(); ++i) {
                                     for (unsigned int j = i + 1; j < p4_2.size(); ++j) {
                                         if (p4_1[i].pt() < 0.0 || p4_2[j].pt() < 0.0)
                                             continue; 
                                         /// need opposite sign dimuons
                                         if ( particle_charges[goodmuons_index[i]] + particle_charges[goodmuons_index[j]] != 0 ) {
                                             continue;
                                         }
                                         /// Add dimuon mass window, 70,110
                                         if ( (p4_1[i] + p4_2[j]).mass() < 70 || (p4_1[i] + p4_2[j]).mass() > 110 ) {
                                             continue;
                                         }
                                         if ( p4_1[i].pt() + p4_2[j].pt() > ptsum) {
                                             ptsum = p4_1[i].pt() + p4_2[j].pt();
                                             if ( p4_1[i].pt() > p4_1[j].pt() ) {
                                                index1 = goodmuons_index[i];
                                                index2 = goodmuons_index[j];
                                             } else {
                                                index1 = goodmuons_index[j];
                                                index2 = goodmuons_index[i]; /// need to return the index1 and index2 as goodmuons pair collection.
                                             }
                                         }
                                     }
                                 }
                                 ROOT::RVec<int> DiMuonPair = {index1, index2};
                                 return DiMuonPair;
                                 ///p4_dimuon = p4_dileptonsystem[0];
                                 ///return p4_dimuon; /// return dimuon_pair_p4 order by pt
                             };
    auto df1 = 
        df.Define(outputname, pair_calc_p4byPt, {particle_pts, particle_etas, particle_phis, particle_masses, particle_charges, goodmuons_index});
    return df1;
}
///
///
/// function  to make a flag that if exist dimuon pair in control region
ROOT::RDF::RNode DiMuonFromCR(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &dimuons_index) {
    auto ZCand_Flag = [](const ROOT::RVec<int> &dimuons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dimuons_index.at(0) == -1 || dimuons_index.at(1) == -1 ) {
                                    return 0;
                                 } else { // else exist a dimuon pair that may from Higgs
                                    return 1;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, ZCand_Flag, {dimuons_index});
    return df1;
}
///
/// function to use dimuon pair index to get dmuon p4 in control region
ROOT::RDF::RNode ZControlDiMuonPairP4(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &particle_pts,
                                 const std::string &particle_etas,
                                 const std::string &particle_phis,
                                 const std::string &particle_masses,
                                 const std::string &dimuons_index) {
    auto dimuon_calc_p4byPt = [](const ROOT::RVec<float> &particle_pts,
                                 const ROOT::RVec<float> &particle_etas,
                                 const ROOT::RVec<float> &particle_phis,
                                 const ROOT::RVec<float> &particle_masses,
                                 const ROOT::RVec<int> &dimuons_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( dimuons_index.at(0) == -1 || dimuons_index.at(1) == -1 ) {
                                    return ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float,default_float);
                                 } else {
                                    // int leading    mu1 = dimuons_index[0];
                                    // int subleading mu2 = dimuons_index[1];
                                    for (unsigned int k = 0; k < (int)dimuons_index.size(); ++k) {
                                        p4.push_back(ROOT::Math::PtEtaPhiMVector(particle_pts.at(dimuons_index[k]),
                                                                                particle_etas.at(dimuons_index[k]),
                                                                                particle_phis.at(dimuons_index[k]),
                                                                                particle_masses.at(dimuons_index[k])));
                                    }
                                    p4_dimuon = ROOT::Math::PtEtaPhiMVector((p4[0]+p4[1]).pt(),
                                                                            (p4[0]+p4[1]).eta(),
                                                                            (p4[0]+p4[1]).phi(),
                                                                            (p4[0]+p4[1]).mass());
                                    return p4_dimuon;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, dimuon_calc_p4byPt, {particle_pts, particle_etas, particle_phis, particle_masses, dimuons_index});
    return df1;
}
///
/// function to pick ele muon pair in Top control region
ROOT::RDF::RNode TOP_EleMuPair_CR(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_pts,
                                 const std::string &muon_etas,
                                 const std::string &muon_phis,
                                 const std::string &muon_masses,
                                 const std::string &muon_charges,
                                 const std::string &goodmuons_index,
                                 const std::string &ele_pts,
                                 const std::string &ele_etas,
                                 const std::string &ele_phis,
                                 const std::string &ele_masses,
                                 const std::string &ele_charges,
                                 const std::string &baseeles_index) {
    auto pair_calc_p4byPt = [](const ROOT::RVec<float> &muon_pts,
                               const ROOT::RVec<float> &muon_etas,
                               const ROOT::RVec<float> &muon_phis,
                               const ROOT::RVec<float> &muon_masses,
                               const ROOT::RVec<int> &muon_charges,
                               const ROOT::RVec<int> &goodmuons_index,
                               const ROOT::RVec<float> &ele_pts,
                               const ROOT::RVec<float> &ele_etas,
                               const ROOT::RVec<float> &ele_phis,
                               const ROOT::RVec<float> &ele_masses,
                               const ROOT::RVec<int> &ele_charges,
                               const ROOT::RVec<int> &baseeles_index) {
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_mu;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4_ele;
                                 for (unsigned int k = 0; k < (int)goodmuons_index.size(); ++k) {
                                    try {
                                        p4_mu.push_back(ROOT::Math::PtEtaPhiMVector(muon_pts.at(goodmuons_index[k]),  ///goodmuons_index[k] points to the good muon index k
                                                                         muon_etas.at(goodmuons_index[k]),          // k = 0, points to goodmuon_index[0]
                                                                         muon_phis.at(goodmuons_index[k]),          // k ,points to goodmuon_index[k]
                                                                         muon_masses.at(goodmuons_index[k])));      // index what I want is goodmuon_index[k] k,i or j
                                    } catch (const std::out_of_range &e) {
                                        p4_mu.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 for (unsigned int k = 0; k < (int)baseeles_index.size(); ++k) {
                                    try {
                                        p4_ele.push_back(ROOT::Math::PtEtaPhiMVector(ele_pts.at(baseeles_index[k]),  ///goodmuons_index[k] points to the good muon index k
                                                                         ele_etas.at(baseeles_index[k]),          // k = 0, points to goodmuon_index[0]
                                                                         ele_phis.at(baseeles_index[k]),          // k ,points to goodmuon_index[k]
                                                                         ele_masses.at(baseeles_index[k])));      // index what I want is goodmuon_index[k] k,i or j
                                    } catch (const std::out_of_range &e) {
                                        p4_ele.push_back(ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float));
                                    }
                                 }
                                 float ptsum = -1;
                                 int index1 = -1,index2 = -1;
                                 for (unsigned int i = 0; i < p4_mu.size(); ++i) {
                                     for (unsigned int j = 0; j < p4_ele.size(); ++j) {
                                         if (p4_mu[i].pt() < 0.0 || p4_ele[j].pt() < 0.0)
                                             continue; 
                                         /// need opposite sign ele and muon
                                         if ( muon_charges[goodmuons_index[i]] + ele_charges[baseeles_index[j]] != 0 ) {
                                             continue;
                                         }
                                         /// Add dimuon mass window, 70,110,150
                                         if ( (p4_mu[i] + p4_ele[j]).mass() < 110 || (p4_mu[i] + p4_ele[j]).mass() > 150 ) {
                                             continue;
                                         }
                                         if ( p4_mu[i].pt() + p4_ele[j].pt() > ptsum) {
                                             ptsum = p4_mu[i].pt() + p4_ele[j].pt();
                                             index1 = goodmuons_index[i];
                                             index2 = baseeles_index[j];
                                         }
                                     }
                                 }
                                 ROOT::RVec<int> EleMuPair = {index1, index2};
                                 return EleMuPair;
                                 ///p4_dimuon = p4_dileptonsystem[0];
                                 ///return p4_dimuon; /// return dimuon_pair_p4 order by pt
                             };
    auto df1 = 
        df.Define(outputname, pair_calc_p4byPt, {muon_pts, muon_etas, muon_phis, muon_masses, muon_charges, goodmuons_index,ele_pts, ele_etas, ele_phis, ele_masses, ele_charges, baseeles_index});
    return df1;
}
///
/// function  to make a flag that if exist ele mu pair in Top control region
ROOT::RDF::RNode EleMuFromCR(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &elemu_index) {
    auto EleMuCand_Flag = [](const ROOT::RVec<int> &elemu_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_dimuon;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( elemu_index.at(0) == -1 || elemu_index.at(1) == -1 ) {
                                    return 0;
                                 } else { // else exist a dimuon pair that may from Higgs
                                    return 1;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, EleMuCand_Flag, {elemu_index});
    return df1;
}
///
/// function to use ele mu pair index to get ele mu p4 in top control region
/// elemu_index, index[0] stands mu, index[1] stands ele
ROOT::RDF::RNode TopControlEleMuPairP4(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &muon_pts,
                                 const std::string &muon_etas,
                                 const std::string &muon_phis,
                                 const std::string &muon_masses,
                                 const std::string &ele_pts,
                                 const std::string &ele_etas,
                                 const std::string &ele_phis,
                                 const std::string &ele_masses,
                                 const std::string &elemu_index) {
    auto elemu_calc_p4byPt = [](const ROOT::RVec<float> &muon_pts,
                                 const ROOT::RVec<float> &muon_etas,
                                 const ROOT::RVec<float> &muon_phis,
                                 const ROOT::RVec<float> &muon_masses,
                                 const ROOT::RVec<float> &ele_pts,
                                 const ROOT::RVec<float> &ele_etas,
                                 const ROOT::RVec<float> &ele_phis,
                                 const ROOT::RVec<float> &ele_masses,
                                 const ROOT::RVec<int> &elemu_index) {
                                 ROOT::Math::PtEtaPhiMVector p4_elemu;
                                 std::vector<ROOT::Math::PtEtaPhiMVector> p4;
                                 if ( elemu_index.at(0) == -1 || elemu_index.at(1) == -1 ) {
                                    return ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float,default_float);
                                 } else {
                                    // int  mu  = elemu_index[0];
                                    // int  ele = elemu_index[1];
                                    p4.push_back(ROOT::Math::PtEtaPhiMVector(muon_pts.at(elemu_index[0]),
                                                                             muon_etas.at(elemu_index[0]),
                                                                             muon_phis.at(elemu_index[0]),
                                                                             muon_masses.at(elemu_index[0])));
                                    p4.push_back(ROOT::Math::PtEtaPhiMVector(ele_pts.at(elemu_index[1]),
                                                                             ele_etas.at(elemu_index[1]),
                                                                             ele_phis.at(elemu_index[1]),
                                                                             ele_masses.at(elemu_index[1])));
                                    p4_elemu = ROOT::Math::PtEtaPhiMVector((p4[0]+p4[1]).pt(),
                                                                            (p4[0]+p4[1]).eta(),
                                                                            (p4[0]+p4[1]).phi(),
                                                                            (p4[0]+p4[1]).mass());
                                    return p4_elemu;
                                 }
                             };
    auto df1 = 
        df.Define(outputname, elemu_calc_p4byPt, {muon_pts, muon_etas, muon_phis, muon_masses, ele_pts, ele_etas, ele_phis, ele_masses, elemu_index});
    return df1;
}
///
/// function  to make a flag that if exist dimuon pair from Higgs
ROOT::RDF::RNode LeadingFatJetSoftDropMass(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &fatjet_softdrop_mass, 
                                 const std::string &good_fatjets_index) {
    auto FatJetSDM = [](const ROOT::RVec<float> &fatjet_softdrop_mass,
                              const ROOT::RVec<int> &good_fatjets_index) {
                                // just return the MassSD after good fatjet
                                 float FatJet_MSD = fatjet_softdrop_mass.at(good_fatjets_index[0]);
                                 return FatJet_MSD;
                             };
    auto df1 = 
        df.Define(outputname, FatJetSDM, {fatjet_softdrop_mass, good_fatjets_index});
    return df1;
}
///
ROOT::RDF::RNode LeadingFatJetVar(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &fatjet_deepTag, 
                                 const std::string &good_fatjets_index) {
    auto FatJetVar = [](const ROOT::RVec<float> &fatjet_deepTag,
                              const ROOT::RVec<int> &good_fatjets_index) {
                                // just return the MassSD after good fatjet
                                 float FatJet_deepTag = fatjet_deepTag.at(good_fatjets_index[0]);
                                 return FatJet_deepTag;
                             };
    auto df1 = 
        df.Define(outputname, FatJetVar, {fatjet_deepTag, good_fatjets_index});
    return df1;
}

///
/// Make Higgs To MuMu Pair Return to a mask
// ROOT::RDF::RNode HiggsToMuMu_Cand(ROOT::RDF::RNode df, const std::string &maskname,
//                         const std::string &dimuon_p4) {
//     // std::string dimuon_pts = dimuon_p4.pt();
//     std::cout << "************************" << std::endl;
//     std::cout << "dimuon_p4" << dimuon_p4 << std::endl;
//     ROOT::Math::PtEtaPhiMVector dimuon_p4_vec;
//     std::istringstream iss(dimuon_p4);
//     iss >> dimuon_p4_vec;
//     std::cout << "dimuon_p4_vec.Pt() is" << dimuon_p4_vec.Pt() << std::endl;
//     // std::cout << "dimuon_p4.at(1) is" << dimuon_p4.Eta() << std::endl;
//     // std::cout << "dimuon_p4.at(2) is" << dimuon_p4.Phi() << std::endl;
//     // std::cout << "dimuon_p4.at(3) is" << dimuon_p4.M() << std::endl;
//     const float dimuon_pt = dimuon_p4_vec.Pt(); // maybe 0123 stands for pt eta phi mass?
//     auto lambda = [dimuon_pt]() {
//         int mask = ( dimuon_pt > 0 ? 1 : 0 ) ;
//         return mask;
//     };
//     auto df1 = df.Define(maskname, lambda, {});
//     return df1;
// }
///
// ROOT::RDF::RNode GetFirstElement(ROOT::RDF::RNode df,
//                                  const std::vector<ROOT::Math::PtEtaPhiMVector> &input_p4_vec,
//                                  const std::string &outputname){
//     auto FirstEle = [](const float &pt,
//                        const float &eta,
//                        const float &phi,
//                        const float &mass) { 
//         ROOT::Math::PtEtaPhiMVector output_p4;
//         output_p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
//         return output_p4; 
//         };
//     return df.Define(outputname, FirstEle, {input_p4_vec[0].pt(), input_p4_vec[0].eta(), input_p4_vec[0].phi(), input_p4_vec[0].mass()});                                
// }
/// end write by botao


/// [begin] a couple of general functions added in developing vhmm analysis
///
/// Function to select objects above a threshold on a variable, using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the varible column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] the lower threshold on the variable
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutVarMin(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMin(threshold), {quantity});
    return df1;
}
/// Function to select objects below a threshold on a variable, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the variable column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] the upper threshold on the variable
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutVarMax(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMax(threshold), {quantity});
    return df1;
}
ROOT::RDF::RNode CutVarMaxUChar(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const unsigned char &threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMaxUChar(threshold), {quantity});
    return df1;
}
/// [end] a couple of general functions added in developing vhmm analysis

/// Function to select objects above a pt threshold, using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the pt column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] ptThreshold minimal pt value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutPt(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &ptThreshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMin(ptThreshold), {quantity});
    return df1;
}
/// Function to select objects blow an eta threshold, using
/// basefunctions::FilterAbsMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the eta column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] EtaThreshold maximal eta value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutEta(ROOT::RDF::RNode df, const std::string &quantity,
                        const std::string &maskname,
                        const float &EtaThreshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterAbsMax(EtaThreshold),
                         {quantity});
    return df1;
}
/// Function to select objects below an Dz threshold, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the Dz column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] Threshold maximal Dz value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDz(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &Threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterAbsMax(Threshold), {quantity});
    return df1;
}
/// Function to select objects below an Dxy threshold, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the Dxy column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] Threshold maximal Dxy value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDxy(ROOT::RDF::RNode df, const std::string &quantity,
                        const std::string &maskname, const float &Threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterAbsMax(Threshold), {quantity});
    return df1;
}
/// Function to select objects with eta dependent upper and lower thesholds
/// for a given variable
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the mask to be added as column to the
/// \param[in] etaColumnName name of the eta column in the NanoAOD dataframe
/// \param[in] cutVarColumnName name of the variable column to apply the
/// selection in the NanoAOD dataframe \param[in] etaBoundary boundary of
/// absolute eta for the barrel and endcap regions of the detector \param[in]
/// lowerThresholdBarrel lower threshold for the barrel \param[in]
/// upperThresholdBarrel upper threshold for the barrel \param[in]
/// lowerThresholdEndcap lower threshold for the endcap \param[in]
/// upperThresholdEndcap upper threshold for the barrel
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutVariableBarrelEndcap(
    ROOT::RDF::RNode df, const std::string &maskname,
    const std::string &etaColumnName, const std::string &cutVarColumnName,
    const float &etaBoundary, const float &lowerThresholdBarrel,
    const float &upperThresholdBarrel, const float &lowerThresholdEndcap,
    const float &upperThresholdEndcap) {
    auto lambda = [etaBoundary, lowerThresholdBarrel, upperThresholdBarrel,
                   lowerThresholdEndcap,
                   upperThresholdEndcap](const ROOT::RVec<float> &eta,
                                         const ROOT::RVec<float> &variable) {
        ROOT::RVec<int> mask =
            (((abs(eta) < etaBoundary) && (variable >= lowerThresholdBarrel) &&
              (variable < upperThresholdBarrel)) ||
             ((abs(eta) >= etaBoundary) && (variable >= lowerThresholdEndcap) &&
              (variable < upperThresholdEndcap)));
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {etaColumnName, cutVarColumnName});
    return df1;
}

/// Function to take a mask and create a new one where a particle candidate is
/// set to false
///
/// \param[in] df the input dataframe
/// \param[out] outputmaskname the name of the new mask to be added as column to
/// the dataframe \param[in] inputmaskname the name of the input mask \param[in]
/// dileptonpair name of the column of the dileptonpair \param[in] index index
/// of the particle candidate to be ignored by mask
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode VetoCandInMask(ROOT::RDF::RNode df,
                                const std::string &outputmaskname,
                                const std::string &inputmaskname,
                                const std::string &dileptonpair,
                                const int index) {
    return df.Define(outputmaskname,
                     [index, inputmaskname](const ROOT::RVec<int> &mask,
                                            const ROOT::RVec<int> &pair) {
                         Logger::get("VetoCandInMask")
                             ->debug("Vetoing the selected candidate (index "
                                     "{}) from the mask {}",
                                     index, inputmaskname);
                         auto newmask = mask;
                         if (pair.at(index) >= 0)
                             newmask.at(pair.at(index)) = 0;
                         return newmask;
                     },
                     {inputmaskname, dileptonpair});
}

/// Function to filter events based on a mask. If the mask contains at least
/// one object, the event is filtered
///
///   \param df the input dataframe
///   \param maskname the name of the column containing the vetomap
///    to be used
///
///   \return a new df with the events filtered
ROOT::RDF::RNode FilterMasks(ROOT::RDF::RNode df, const std::string &maskname) {
    auto df1 = df.Filter(
        [](const ROOT::RVec<Int_t> &mask) {
            auto result = Any(mask);
            return result;
        },
        {maskname});
    return df1;
}

/// Function to generate a veto based on a mask. If the mask contains at least
/// one object, the veto flag is set, meaning that this event contains at least
/// one object matching the requirements of the veto map.
///
/// \code
///  In the veto column
///  1 --> the event contains at least one object matching the veto map
///  0 --> the event does not contain any object matching the veto map
/// \endcode
///
/// \param df the input dataframe
/// \param outputname name of the new quantity containing the veto flags
/// \param vetomap the name of the column containing the vetomap to be used
///
/// \return a new df containing the veto flag column
ROOT::RDF::RNode LeptonVetoFlag(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &vetomap) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() != 0;
                     },
                     {vetomap});
}

/// Function to create a boolian flag based on the number of non-zero masks.
/// The flag is set to true if the number of non-zero masks is zero.
///
/// \param df the input dataframe
/// \param outputname the name of the output column that is created
/// \param vetomap the name of the column containing the input mask to be used
///
/// \return a new df containing the output flag column
ROOT::RDF::RNode IsEmptyFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &vetomap) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() == 0;
                     },
                     {vetomap});
}

/// Function to create a boolian flag based on the number of non-zero masks.
/// The flag is set to true if the number of non-zero masks matches the
/// allowed value provided as an input.
///
/// \param df the input dataframe
/// \param outputname the name of the output column that is created
/// \param map the name of the column containing the input mask to be used
/// \param n the allowed number of non-zero masks
///
/// \return a new df containing the output flag column
ROOT::RDF::RNode CutNFlag(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &map, const int &n) {
    return df.Define(outputname,
                     [n](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() == n;
                     },
                     {map});
}

/// Function to create a column for a vector of indices of objects based on
/// input masks. Indices of objects with non-zero mask are stored in the output
/// column.
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] inputmaskname the name of input masks
///
/// \return a dataframe containing a vector of indices for the selected objects
ROOT::RDF::RNode SelectedObjects(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &inputmaskname) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         Logger::get("SelectedObjects")
                             ->debug("size = {}",
                                     ROOT::VecOps::Nonzero(mask).size());
                         return static_cast<ROOT::VecOps::RVec<int>>(
                             ROOT::VecOps::Nonzero(mask));
                     },
                     {inputmaskname});
}

/**
 * @brief Function used to veto a particle, if it is overlapping within a given
 * DeltaR value with another particle within a mask. For the particle to test,
 * the lorentz vector is used, for the veto, the mask and all input quantities
 * needed to calculate the lorentz vector are used. Returns true if the particle
 * is vetoed, false if it is not.
 *
 * @param df The input dataframe
 * @param output_flag The name of the veto flag to be added as column to the
 * dataframe
 * @param p4 The name of the Lorentz vector column to be used for the particle
 * to test
 * @param particle_mask The name of the mask to be used for the veto
 * @param particle_pt The name of the pt column to be used for the particles to
 * test against
 * @param particle_eta The name of the eta column to be used for the particles
 * to test against
 * @param particle_phi The name of the phi column to be used for the particles
 * to test against
 * @param particle_mass The name of the mass column to be used for the particles
 * to test against
 * @param dR_cut The maximum dR to be used for the veto
 * @return ROOT::RDF::RNode a dataframe containing the new veto flag
 */

ROOT::RDF::RNode DeltaRParticleVeto(
    ROOT::RDF::RNode df, const std::string &output_flag, const std::string &p4,
    const std::string &particle_mask, const std::string &particle_pt,
    const std::string &particle_eta, const std::string &particle_phi,
    const std::string &particle_mass, const float dR_cut) {
    auto veto_overlapping_particle =
        [dR_cut](const ROOT::Math::PtEtaPhiMVector &p4,
                 const ROOT::RVec<float> &particle_pt,
                 const ROOT::RVec<float> &particle_eta,
                 const ROOT::RVec<float> &particle_phi,
                 const ROOT::RVec<float> &particle_mass,
                 const ROOT::RVec<int> &particle_mask) {
            // for all particles in the mask, check if they overlap with the
            // particle, if so, return true
            ROOT::RVec<int> valid_particle_indices =
                ROOT::VecOps::Nonzero(particle_mask);
            const auto selected_pt =
                ROOT::VecOps::Take(particle_pt, valid_particle_indices);
            const auto selected_eta =
                ROOT::VecOps::Take(particle_eta, valid_particle_indices);
            const auto selected_phi =
                ROOT::VecOps::Take(particle_phi, valid_particle_indices);
            const auto selected_mass =
                ROOT::VecOps::Take(particle_mass, valid_particle_indices);
            auto selected_p4s =
                ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
                    selected_pt, selected_eta, selected_phi, selected_mass);
            for (const auto &p4_test : selected_p4s) {
                if (ROOT::Math::VectorUtil::DeltaR(p4_test, p4) < dR_cut) {
                    return true;
                }
            }
            // if no particle is close enough to the p4, return false
            return false;
        };
    auto df1 = df.Define(output_flag, veto_overlapping_particle,
                         {p4, particle_pt, particle_eta, particle_phi,
                          particle_mass, particle_mask});
    return df1;
}

/// Function to correct object mass in alignment with object pt correction
///
/// \param[in] df the input dataframe
/// \param[out] corrected_mass the name of the corrected masses to be determined
/// \param[in] raw_mass name of the input mass \param[in] raw_pt name of the
/// uncorrected object pts \param[in] corrected_pt name of the corrected object
/// pts
///
/// \return a dataframe containing the modified object masses
ROOT::RDF::RNode ObjectMassCorrectionWithPt(ROOT::RDF::RNode df,
                                            const std::string &corrected_mass,
                                            const std::string &raw_mass,
                                            const std::string &raw_pt,
                                            const std::string &corrected_pt) {
    auto mass_correction_lambda =
        [](const ROOT::RVec<float> &mass_values,
           const ROOT::RVec<float> &pt_values,
           const ROOT::RVec<float> &corrected_pt_values) {
            ROOT::RVec<float> corrected_mass_values(mass_values.size());
            for (int i = 0; i < mass_values.size(); i++) {
                corrected_mass_values[i] = mass_values.at(i) *
                                           corrected_pt_values.at(i) /
                                           pt_values.at(i);
            }
            return corrected_mass_values;
        };
    auto df1 = df.Define(corrected_mass, mass_correction_lambda,
                         {raw_mass, raw_pt, corrected_pt});
    return df1;
}

/// Function to check whether at least one lepton pair is present
///
/// \param[in] df the input dataframe
/// \param[out] output_flag the name of the bool column that is created
/// \param[in] leptons_pt name of the input pt column of the lepton collection
/// \param[in] leptons_eta name of the input eta column of the lepton collection
/// \param[in] leptons_phi name of the input phi column of the lepton collection
/// \param[in] leptons_mass name of the input mass column of the lepton
/// collection \param[in] leptons_charge name of the input charge column of the
/// lepton collection \param[in] leptons_mask name of the input mask column of
/// the lepton collection that marks lepton to be taken into account \param[in]
/// dR_cut minimum required angular distance between the leptons
///
/// \return a dataframe containing the new bool column
ROOT::RDF::RNode CheckForDiLeptonPairs(
    ROOT::RDF::RNode df, const std::string &output_flag,
    const std::string &leptons_pt, const std::string &leptons_eta,
    const std::string &leptons_phi, const std::string &leptons_mass,
    const std::string &leptons_charge, const std::string &leptons_mask,
    const float dR_cut) {
    auto pair_finder_lambda = [dR_cut](const ROOT::RVec<float> &pt_values,
                                       const ROOT::RVec<float> &eta_values,
                                       const ROOT::RVec<float> &phi_values,
                                       const ROOT::RVec<float> &mass_values,
                                       const ROOT::RVec<int> &charge_values,
                                       const ROOT::RVec<int> &mask) {
        const auto valid_lepton_indices = ROOT::VecOps::Nonzero(mask);
        for (auto it1 = valid_lepton_indices.begin();
             it1 != valid_lepton_indices.end(); it1++) {
            for (auto it2 = it1 + 1; it2 != valid_lepton_indices.end(); it2++) {
                if (charge_values.at(*it1) != charge_values.at(*it2)) {
                    auto p4_1 = ROOT::Math::PtEtaPhiMVector(
                        pt_values.at(*it1), eta_values.at(*it1),
                        phi_values.at(*it1), mass_values.at(*it1));
                    auto p4_2 = ROOT::Math::PtEtaPhiMVector(
                        pt_values.at(*it2), eta_values.at(*it2),
                        phi_values.at(*it2), mass_values.at(*it2));
                    if (ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2) >= dR_cut)
                        return true;
                }
            }
        }
        return false;
    };
    auto df1 = df.Define(output_flag, pair_finder_lambda,
                         {leptons_pt, leptons_eta, leptons_phi, leptons_mass,
                          leptons_charge, leptons_mask});
    return df1;
}

/// Muon specific functions
namespace muon {
/// Function to cut on muons based on the muon ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to cut on muons based on the muon isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}
/// Function to cut on muons based on the muon isolation using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold minimal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode AntiCutIsolation(ROOT::RDF::RNode df,
                                  const std::string &maskname,
                                  const std::string &isolationName,
                                  const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMin(Threshold),
                         {isolationName});
    return df1;
}
/// Function to cut on muons based on the muon signature: isTracker or isGlobal
///
/// \param[in] df the input dataframe
/// \param[in] isTracker name of the signature column in the NanoAOD
/// \param[in] isGlobal name of the signature column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe 
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsTrackerOrIsGlobal(ROOT::RDF::RNode df, const std::string &isTracker, const std::string &isGlobal, const std::string &maskname) {
        auto lambda = [](const ROOT::RVec<Bool_t>  &tracker, const ROOT::RVec<Bool_t>  &global) {
            ROOT::RVec<int> mask = (tracker == 1 || global == 1);
            Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                ->debug("istracker {}, isglobal {}, mask {}", tracker, global, mask);
            return mask;
        };
        auto df1 = df.Define(maskname, lambda, {isTracker, isGlobal});
        return df1;
}
/// Function to create a column of vector of random numbers between 0 and 1
/// with size of the input object collection
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] objCollection the name of the input object collection
/// \param[in] seed the seed of the random number generator
///
/// \return a dataframe with the new column
ROOT::RDF::RNode GenerateRndmRVec(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &objCollection, int seed) {
    gRandom->SetSeed(seed);
    auto lambda = [](const ROOT::RVec<int> &objects) {
        const int len = objects.size();
        float rndm[len];
        gRandom->RndmArray(len, rndm);
        ROOT::RVec<float> out = {};
        for (auto &x : rndm) {
            out.push_back(x);
        }
        return out;
    };
    return df.Define(outputname, lambda, {objCollection});
}

/// Function to create a column of Rochester correction applied transverse
/// momentum for data, see https://gitlab.cern.ch/akhukhun/roccor
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] filename the name of Rochester correction file
/// \param[in] position index of the position in the input vector
/// \param[in] objCollection the name of the column containing the input vector
/// \param[in] chargColumn the name of the column containing muon charges
/// \param[in] ptColumn the name of the column containing muon charge values
/// \param[in] etaColumn the name of the column containing muon eta values
/// \param[in] phiColumn the name of the column containing muon phi values
/// \param[in] error_set the error set number
/// \param[in] error_member the error member number
///
/// \return a dataframe with the new column
ROOT::RDF::RNode
applyRoccoRData(ROOT::RDF::RNode df, const std::string &outputname,
                const std::string &filename, const int &position,
                const std::string &objCollection,
                const std::string &chargColumn, const std::string &ptColumn,
                const std::string &etaColumn, const std::string &phiColumn,
                int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set,
                   error_member](const ROOT::RVec<int> &objects,
                                 const ROOT::RVec<int> &chargCol,
                                 const ROOT::RVec<float> &ptCol,
                                 const ROOT::RVec<float> &etaCol,
                                 const ROOT::RVec<float> &phiCol) {
        const int index = objects.at(position);
        double pt_rc =
            ptCol.at(index) * rc.kScaleDT(chargCol.at(index), ptCol.at(index),
                                          etaCol.at(index), phiCol.at(index),
                                          error_set, error_member);
        return pt_rc;
    };

    return df.Define(
        outputname, lambda,
        {objCollection, chargColumn, ptColumn, etaColumn, phiColumn});
}

/// Function to create a column of Rochester correction applied transverse
/// momentum for MC, see https://gitlab.cern.ch/akhukhun/roccor
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] filename the name of Rochester correction file
/// \param[in] position index of the position in the input vector
/// \param[in] objCollection the name of the column containing the input vector
/// \param[in] chargColumn the name of the column containing muon charges
/// \param[in] ptColumn the name of the column containing muon charge values
/// \param[in] etaColumn the name of the column containing muon eta values
/// \param[in] phiColumn the name of the column containing muon phi values
/// \param[in] genPtColumn the name of the column containing gen-level
/// transverse momentum value of the target muon \param[in] nTrackerLayersColumn
/// the name of the column containing number of tracker layers values \param[in]
/// rndmColumn the name of the column containing random number generated for
/// each muon \param[in] error_set the error set number \param[in] error_member
/// the error member number
///
/// \return a dataframe with the new column
ROOT::RDF::RNode
applyRoccoRMC(ROOT::RDF::RNode df, const std::string &outputname,
              const std::string &filename, const int &position,
              const std::string &objCollection, const std::string &chargColumn,
              const std::string &ptColumn, const std::string &etaColumn,
              const std::string &phiColumn, const std::string &genPtColumn,
              const std::string &nTrackerLayersColumn,
              const std::string &rndmColumn, int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set, error_member](
                      const ROOT::RVec<int> &objects,
                      const ROOT::RVec<int> &chargCol,
                      const ROOT::RVec<float> &ptCol,
                      const ROOT::RVec<float> &etaCol,
                      const ROOT::RVec<float> &phiCol, const float &genPt,
                      const ROOT::RVec<int> &nTrackerLayersCol,
                      const ROOT::RVec<float> &rndmCol) {
        double pt_rc = default_float;
        const int index = objects.at(position);
        if (genPt > 0.) {
            pt_rc = ptCol.at(index) *
                    rc.kSpreadMC(chargCol.at(index), ptCol.at(index),
                                 etaCol.at(index), phiCol.at(index), genPt,
                                 error_set, error_member);
        } else {
            pt_rc = ptCol.at(index) *
                    rc.kSmearMC(chargCol.at(index), ptCol.at(index),
                                etaCol.at(index), phiCol.at(index),
                                nTrackerLayersCol.at(index),
                                rndmCol.at(position), error_set, error_member);
        }

        return pt_rc;
    };

    return df.Define(outputname, lambda,
                     {objCollection, chargColumn, ptColumn, etaColumn,
                      phiColumn, genPtColumn, nTrackerLayersColumn,
                      rndmColumn});
}
} // end namespace muon
/// Tau specific functions
namespace tau {
/// Function to cut on taus based on the tau decay mode
///
/// \param[in] df the input dataframe
/// \param[in] tau_dms name of the column with tau decay modes
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] SelectedDecayModes a `std::vector<int>` containing the
/// decay modes, that should pass the cut
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDecayModes(ROOT::RDF::RNode df, const std::string &maskname,
                               const std::string &tau_dms,
                               const std::vector<int> &SelectedDecayModes) {
    auto df1 = df.Define(
        maskname,
        [SelectedDecayModes](const ROOT::RVec<Int_t> &decaymodes) {
            ROOT::RVec<int> mask;
            for (auto n : decaymodes) {
                mask.push_back(int(std::find(SelectedDecayModes.begin(),
                                             SelectedDecayModes.end(),
                                             n) != SelectedDecayModes.end()));
            }
            return mask;
        },
        {tau_dms});
    return df1;
}
/// Function to cut taus based on the tau ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutTauID(ROOT::RDF::RNode df, const std::string &maskname,
                          const std::string &nameID, const int &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterID(idxID), {nameID});
    return df1;
}
/// Function to correct e to tau fake pt
///
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] df the input dataframe
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_dm0_b scale factor to be applied to taus with decay mode 0 and
/// eta region barrel \param[in] sf_dm1_b scale factor to be applied to taus
/// with decay mode 1 and eta region barrel \param[in] sf_dm0_e scale factor to
/// be applied to taus with decay mode 0 and eta region endcap \param[in]
/// sf_dm1_e scale factor to be applied to taus with decay mode 1 and eta region
/// endcap name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const std::string &decayMode, const std::string &genMatch,
                     const std::string &sf_file, const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, sf_dm0_b, sf_dm1_b,
                                     sf_dm0_e, sf_dm1_e](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 1 || genmatch.at(i) == 3) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0 &&
                    std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 0 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("ptcorrection ele fake")
                ->debug("tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct mu to tau fake pt
///
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] df the input dataframe
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_es scale factor to be applied to taus faked by muons
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda =
        [evaluator, idAlgorithm, sf_es](const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    // only considering wanted tau decay modes
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_es});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    Logger::get("mu fake")->debug(
                        "tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
                }
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct tau pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt \param[in] decayMode
/// \param[in] sf_dm0 scale factor to be applied to taus with decay mode 0
/// \param[in] sf_dm1 scale factor to be applied to other 1 prong taus
/// \param[in] sf_dm10 scale factor to be applied to taus with decay mode 10
/// \param[in] sf_dm11 scale factor to be applied to other 3 prong taus
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &decayMode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11) {
    auto tau_pt_correction_lambda =
        [sf_dm0, sf_dm1, sf_dm10, sf_dm11](const ROOT::RVec<float> &pt_values,
                                           const ROOT::RVec<int> &decay_modes) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (decay_modes.at(i) == 0)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm0;
                else if (decay_modes.at(i) > 0 && decay_modes.at(i) < 5)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm1;
                else if (decay_modes.at(i) == 10)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm10;
                else if (decay_modes.at(i) > 10 && decay_modes.at(i) < 15)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm11;
                else
                    corrected_pt_values[i] = pt_values.at(i);
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, tau_pt_correction_lambda, {pt, decayMode});
    return df1;
}
/// Function to correct tau pt with correctionlib
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] DM0 variation decay mode 0
/// \param[in] DM1 variation decay mode 1
/// \param[in] DM10 variation decay mode 10
/// \param[in] DM11 variation decay mode 11
/// DM values: "nom","up","down"
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, DM0, DM1, DM10,
                                     DM11](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 5) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM0});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM1});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 10) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM10});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 11) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM11});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("tauEnergyCorrection")
                ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                        pt_values.at(i), corrected_pt_values.at(i),
                        decay_modes.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
} // end namespace tau

namespace electron {
/// Function to correct electron pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected electron pt to be calculated
/// \param[in] pt name of the raw electron pt
/// \param[in] eta the name of the raw electron eta
/// \param[in] sf_barrel scale factor to be applied to electrons in the barrel
/// \param[in] sf_endcap scale factor to be applied to electrons in the endcap
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap) {
    auto electron_pt_correction_lambda =
        [sf_barrel, sf_endcap](const ROOT::RVec<float> &pt_values,
                               const ROOT::RVec<float> &eta) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (abs(eta.at(i)) <= 1.479) {
                    corrected_pt_values[i] = pt_values.at(i) * sf_barrel;
                } else if (abs(eta.at(i)) > 1.479) {
                    corrected_pt_values[i] = pt_values.at(i) * sf_endcap;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, electron_pt_correction_lambda, {pt, eta});
    return df1;
}
/// Function to cut electrons based on the electron MVA ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
ROOT::RDF::RNode CutUCharID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID, const unsigned char &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterJetUCharID(idxID), {nameID});
    return df1;
}
/// Function to cut electrons based on the cut based electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe \param[in] nameID name of the ID column in the NanoAOD
/// \param[in] IDvalue value of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const int &IDvalue) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMinInt(IDvalue), {nameID});
    return df1;
}
/// Function to cut electrons based on failing the cut based electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe \param[in] nameID name of the ID column in the NanoAOD
/// \param[in] IDvalue value of the WP the has to be failed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode AntiCutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                             const std::string &nameID, const int &IDvalue) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMaxInt(IDvalue), {nameID});
    return df1;
}

/// Function to cut electrons based on the electron isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}
/// Function to select electrons below an Dxy and Dz threshold, based on the
/// electrons supercluster
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the electron eta column in the NanoAOD
/// \param[in] quantity name of the electron deltaEtaSC column in the NanoAOD
/// \param[in] quantity name of the Dxy column in the NanoAOD
/// \param[in] quantity name of the Dz column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe
/// \param[in] abs(eta) of the EB-EE transition
/// \param[in] Threshold maximal Dxy value in the barrel
/// \param[in] Threshold maximal Dz value in the barrel
/// \param[in] Threshold maximal Dxy value in the endcap
/// \param[in] Threshold maximal Dz value in the endcap
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIP(ROOT::RDF::RNode df, const std::string &eta,
                       const std::string &detasc, const std::string &dxy,
                       const std::string &dz, const std::string &maskname,
                       const float &abseta_eb_ee, const float &max_dxy_eb,
                       const float &max_dz_eb, const float &max_dxy_ee,
                       const float &max_dz_ee) {
    auto lambda = [abseta_eb_ee, max_dxy_eb, max_dz_eb, max_dxy_ee,
                   max_dz_ee](const ROOT::RVec<float> &eta,
                              const ROOT::RVec<float> &detasc,
                              const ROOT::RVec<float> &dxy,
                              const ROOT::RVec<float> &dz) {
        ROOT::RVec<int> mask = (((abs(eta + detasc) < abseta_eb_ee) &&
                                 (dxy < max_dxy_eb) && (dz < max_dz_eb)) ||
                                ((abs(eta + detasc) >= abseta_eb_ee) &&
                                 (dxy < max_dxy_ee) && (dz < max_dz_ee)));
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {eta, detasc, dxy, dz});
    return df1;
}

/// Function to veto electrons in the transition region of EB and EE, based on
/// the electrons supercluster
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the electron eta column in the NanoAOD
/// \param[in] quantity name of the electron deltaEtaSC column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe
/// \param[in] abs(eta) of the beginning of the transition region
/// \param[in] abs(eta) of the end of the transition region
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutGap(ROOT::RDF::RNode df, const std::string &eta,
                        const std::string &detasc, const std::string &maskname,
                        const float &end_eb, const float &start_ee) {
    auto lambda = [end_eb, start_ee](const ROOT::RVec<float> &eta,
                                     const ROOT::RVec<float> &detasc) {
        ROOT::RVec<int> mask =
            (abs(eta + detasc) < end_eb) || (abs(eta + detasc) > start_ee);
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {eta, detasc});
    return df1;
}

} // end namespace electron

} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */
