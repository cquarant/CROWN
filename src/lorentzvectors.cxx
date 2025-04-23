#ifndef GUARDLVECS_H
#define GUARDLVECS_H

#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>

/// Namespace used for lorentzvector operations

namespace lorentzvectors {


ROOT::RDF::RNode buildSFMass(
    ROOT::RDF::RNode df, 
    const std::vector<std::string> &obj_quantities,
    const int pairindex, 
    const std::string &obj_sfmass_name
    ){
        auto df1 = df.Define(
            obj_sfmass_name,
            [pairindex](
                const ROOT::RVec<int> &pair,
                const ROOT::RVec<float> &masses
            ){
                const int index = pair.at(pairindex);
                return masses.at(index);
            }, obj_quantities
        );
        return df1;
    }

ROOT::RDF::RNode build_subjet_postion(
    ROOT::RDF::RNode df, 
    const std::vector<std::string> &obj_quantities,
    const int pairindex, 
    const std::string &obj_sfmass_name
    ){
        auto df1 = df.Define(
            obj_sfmass_name,
            [pairindex](
                const ROOT::RVec<int> &pair,
                const ROOT::RVec<short> &inn
            ){
                const int index = pair.at(pairindex);
                return inn.at(index);
            }, obj_quantities
        );
        return df1;
    }

ROOT::RDF::RNode matchSubJet(
    ROOT::RDF::RNode df, 
    const std::vector<std::string> &obj_quantities,
    const int pairindex0, const int pairindex1, 
    const std::string &obj_sfmass_name
    ){
        auto df1 =df.Define(
            obj_sfmass_name,
            [pairindex0, pairindex1](
                const ROOT::RVec<int> &pair,
                const short &index_from_Fatjet0,
                const short &index_from_Fatjet1
            ){
                int matching = 0;

                const int index0 = pair.at(pairindex0);
                const int index1 = pair.at(pairindex1);
                int index00 = index_from_Fatjet0;
                int index11 = index_from_Fatjet1;
                if (index0==index00 && index1==index11){
                    matching = 1;
                }
                else if (index0==index11 && index1==index00){
                    matching = 1;
                }
                return matching;
            }, obj_quantities
        );
        return df1;
    }

ROOT::RDF::RNode buildSafe999(ROOT::RDF::RNode df, const std::string &outputname) {
    auto df1 = df.Define(
        outputname,
        []() {
            constexpr float safe_value = -999.0;
            Logger::get("lorentzvectors")->debug("Returning safe default mass: {}", safe_value);
            return safe_value;
        }
    );
    return df1;
}

ROOT::RDF::RNode buildSafem1(ROOT::RDF::RNode df, const std::string &outputname) {
    auto df1 = df.Define(
        outputname,
        []() {
            const int safe_value = -1;
            return safe_value;
        }
    );
    return df1;
}
    

/// Function to build the lorentzvector from the pt, eta, phi and mass of a
/// particle. This utilizes the [PtEtaPhiMVector from
/// ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31)
///
/// \param df The input dataframe
/// \param quantities This vector contains the names of the columns:
///    -# pair - In this vector the particle index is stored. This index is used
///    to select the correct particle from the four particle quantity vectors.
///    -# pts - In this vector, the pts of the particle are stored.
///    -# etas - In this vector, the etas of the particle are stored.
///    -# phis - In this vector, the phis of the particle are stored.
///    -# masses - In this vector, the masses of the particle are stored.
///      This order has to be kept!
/// \param outputname The name of the output column in the new dataframe
/// \param position The position in the pair vector, which is used to store the
/// index of the particle in the particle quantity vectors.
///
/// \returns a new dataframe, which contains the new lorentz vector
ROOT::RDF::RNode buildparticle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &quantities,
                               const std::string &outputname,
                               const int &position) {
    auto df1 = df.Define(
        outputname,
        [position, outputname](
            const ROOT::RVec<int> &pair, const ROOT::RVec<float> &pts,
            const ROOT::RVec<float> &etas, const ROOT::RVec<float> &phis,
            const ROOT::RVec<float> &masses) {
            // the index of the particle is stored in the pair vector
            ROOT::Math::PtEtaPhiMVector p4;
            Logger::get("lorentzvectors")
                ->debug("starting to build 4vector {}!", outputname);
            try {
                const int index = pair.at(position);
                Logger::get("lorentzvectors")->debug("pair {}", pair);
                Logger::get("lorentzvectors")->debug("pts {}", pts);
                Logger::get("lorentzvectors")->debug("etas {}", etas);
                Logger::get("lorentzvectors")->debug("phis {}", phis);
                Logger::get("lorentzvectors")->debug("masses {}", masses);
                Logger::get("lorentzvectors")->debug("Index {}", index);

                p4 = ROOT::Math::PtEtaPhiMVector(pts.at(index), etas.at(index),
                                                 phis.at(index),
                                                 masses.at(index));
            } catch (const std::out_of_range &e) {
                p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,
                                                 default_float, default_float);
                Logger::get("lorentzvectors")
                    ->debug("Index not found, retuning dummy vector !");
            }
            Logger::get("lorentzvectors")
                ->debug("P4 - Particle {} : {}", position, p4);
            return p4;
        },
        quantities);
    return df1;
}

ROOT::RDF::RNode build_input_index(ROOT::RDF::RNode df,
                               const std::string &index, const std::string &pt,
                               const std::string &eta , const std::string &phi , const std::string &mass ,
                               const std::string &outputname) {
    auto df1 = df.Define(
        outputname,
        [index, outputname](
            const ROOT::RVec<float> &pts,
            const ROOT::RVec<float> &etas, 
            const ROOT::RVec<float> &phis,
            const ROOT::RVec<float> &masses) {
            // the index of the particle is stored in the pair vector
            ROOT::Math::PtEtaPhiMVector p4;
            Logger::get("lorentzvectors")
                ->debug("starting to build 4vector {}!", outputname);
            try {
                const short index = index;
                Logger::get("lorentzvectors")->debug("pts {}", pts);
                Logger::get("lorentzvectors")->debug("etas {}", etas);
                Logger::get("lorentzvectors")->debug("phis {}", phis);
                Logger::get("lorentzvectors")->debug("masses {}", masses);
                Logger::get("lorentzvectors")->debug("Index {}", index);

                p4 = ROOT::Math::PtEtaPhiMVector(pts.at(index), etas.at(index),
                                                 phis.at(index),
                                                 masses.at(index));
            } catch (const std::out_of_range &e) {
                p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,
                                                 default_float, default_float);
                Logger::get("lorentzvectors")
                    ->debug("Index not found, retuning dummy vector !");
            }
            Logger::get("lorentzvectors")
                ->debug("P4 - Particle {} : {}", index, p4);
            return p4;
        },
        {pt, eta, phi, mass});
    return df1;
}

ROOT::RDF::RNode buildSafeP4(ROOT::RDF::RNode df, const std::string &outputname) {
    auto df1 = df.Define(
        outputname,
        []() {
            constexpr float safe_value = 999.0;
            ROOT::Math::PtEtaPhiMVector p4(safe_value, safe_value, safe_value, safe_value);
            Logger::get("lorentzvectors")->debug("Returning safe default 4-vector: {}", p4);
            return p4;
        }
    );
    return df1;
}

/**
 * @brief Function used to construct a 4-vector for a pair particle.
 *
 * @param df the input dataframe
 * @param obj_quantities a vector of strings containing the names of the
 * quantities used for the 4 vector building. As this function is a wrapper for
 * lorentzvectors::buildparticle, the order of the inputs has to match the order
 * used in that function
 * @param pairindex the index of the particle in the pair to be build
 * @param obj_p4_name name of the column containing the 4-vector
 * @return a new df, containing the thw column
 */

ROOT::RDF::RNode build(ROOT::RDF::RNode df,
                       const std::vector<std::string> &obj_quantities,
                       const int pairindex, const std::string &obj_p4_name) {
    Logger::get("lorentzvectors")->debug("Building {}", obj_p4_name);
    for (auto i : obj_quantities)
        Logger::get("lorentzvectors")->debug("Used object quantities {}", i);
    return lorentzvectors::buildparticle(df, obj_quantities, obj_p4_name,
                                         pairindex);
}

/**
 * @brief Function used to construct the missing transverse energy lorentz
 * vector.
 *
 * @param df the input dataframe
 * @param met_pt the met pt value
 * @param met_phi the met phi value
 * @param outputname name of the new column containing the missing transverse
 * energy lorentz vector.
 * @return a new df, containing the new column
 */

ROOT::RDF::RNode buildMet(ROOT::RDF::RNode df, const std::string &met_pt,
                          const std::string &met_phi,
                          const std::string &outputname) {
    auto construct_metvector = [](const float &pt, const float &phi) {
        // for Met, eta is zero
        auto met = ROOT::Math::PtEtaPhiEVector(pt, 0, phi, pt);
        // cast Met vector to a ROOT::Math::PtEtaPhiMVector to make latter
        // functions easier to use
        return (ROOT::Math::PtEtaPhiMVector)met;
    };
    return df.Define(outputname, construct_metvector, {met_pt, met_phi});
}
/**
 * @brief Function to scale a lorentz vector by a scale factor.
 *
 * @param df the input dataframe
 * @param outputname name of the new column containing the missing transverse
 * energy lorentz vector.
 * @param inputvector a vector of the two names of the columns containing the
 * required lorentz vectors
 * @param p4_sf the scale factor, that is applied to the input lorentz vector
 * @return a new df, containing the new column
 */
ROOT::RDF::RNode scaleP4(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::vector<std::string> &inputvector,
                         const float &p4_sf) {
    return df.Define(
        outputname,
        [p4_sf](const ROOT::Math::PtEtaPhiMVector &p4) {
            if (p4.pt() < 0.0)
                return ROOT::Math::PtEtaPhiMVector(
                    default_float, default_float, default_float, default_float);
            auto p4_scaled = p4_sf * p4;
            return p4_scaled;
        },
        inputvector);
}

ROOT::RDF::RNode buildHiggs(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &Genparticle_pdgId,
                                 const std::string &Genparticle_statusFlags,
                                 const std::string &Genparticle_pt,
                                 const std::string &Genparticle_eta,
                                 const std::string &Genparticle_phi,
                                 const std::string &Genparticle_mass
                                 ) {
    auto Higgsp4 = [](const ROOT::RVec<int> &genparticle_pdgId,
        const ROOT::RVec<unsigned short> &genparticle_statusFlags,const ROOT::RVec<float> &genparticle_pt,
        const ROOT::RVec<float> &genparticle_eta,const ROOT::RVec<float> &genparticle_phi,
        const ROOT::RVec<float> &genparticle_mass) {
        ROOT::Math::PtEtaPhiMVector p4;
        p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,default_float, default_float);
        for (unsigned int i = 0; i < (int)genparticle_pdgId.size(); ++i ) {
            // check if the gen particle is H
            // assign pt, eta, phi of higgs to p4 once found
            if ( ( abs(genparticle_pdgId.at(i)) == 25 ) && ( (genparticle_statusFlags.at(i) >> 13) & 1 == 1 ) ) {
                p4 = ROOT::Math::PtEtaPhiMVector(genparticle_pt.at(i),genparticle_eta.at(i),
                genparticle_phi.at(i), genparticle_mass.at(i));
                return p4; 
                }
            }
        return p4; // if not found return the default p4
    };
                    
    auto df1 = 
        df.Define(outputname, Higgsp4, 
            {Genparticle_pdgId, Genparticle_statusFlags, Genparticle_pt,Genparticle_eta,
            Genparticle_phi,Genparticle_mass});
    return df1;
}


} // namespace lorentzvectors
#endif /* GUARDLVECS_H */