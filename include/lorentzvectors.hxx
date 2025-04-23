#ifndef GUARDLVECS_H
#define GUARDLVECS_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/utility.hxx"
#include <Math/Vector4D.h>

namespace lorentzvectors {

ROOT::RDF::RNode buildSFMass(
    ROOT::RDF::RNode df, 
    const std::vector<std::string> &obj_quantities,
    const int pairindex, 
    const std::string &obj_sfmass_name
    );

ROOT::RDF::RNode build_subjet_postion(
    ROOT::RDF::RNode df, 
    const std::vector<std::string> &obj_quantities,
    const int pairindex, 
    const std::string &obj_sfmass_name
    );

ROOT::RDF::RNode matchSubJet(
    ROOT::RDF::RNode df, 
    const std::vector<std::string> &obj_quantities,
    const int pairindex0, const int pairindex1, 
    const std::string &obj_sfmass_name
    );

ROOT::RDF::RNode buildSafe999(ROOT::RDF::RNode df, const std::string &outputname);
ROOT::RDF::RNode buildSafem1(ROOT::RDF::RNode df, const std::string &outputname);

ROOT::RDF::RNode buildparticle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &quantities,
                               const std::string &outputname,
                               const int &position);

ROOT::RDF::RNode build_input_index(ROOT::RDF::RNode df,
                               const std::string &index, const std::string &pt,
                               const std::string &eta , const std::string &phi , const std::string &mass ,
                               const std::string &outputname);
                               
ROOT::RDF::RNode buildSafeP4(ROOT::RDF::RNode df, const std::string &outputname);

ROOT::RDF::RNode build(ROOT::RDF::RNode df,
                       const std::vector<std::string> &obj_quantities,
                       const int pairindex, const std::string &obj_p4_name);
ROOT::RDF::RNode buildMet(ROOT::RDF::RNode df, const std::string &met_pt,
                          const std::string &met_phi,
                          const std::string &outputname);
ROOT::RDF::RNode buildHiggs(ROOT::RDF::RNode df, const std::string &outputname,
                                 const std::string &GenPart_pdgId,
                                 const std::string &GenPart_statusFlags,
                                 const std::string &GenParticle_pt,
                                 const std::string &GenParticle_eta,
                                 const std::string &GenParticle_phi,
                                 const std::string &GenParticle_mass
                                 );
/// This function constructs a vectorial sum of n lorentz vectors. Since it uses
/// an aribtrary number of input vectors, it is implemented in the header file.
///
/// \param df The input dataframe
/// \param outputflag The name of the new column
/// \param p4s Parameter pack of column names that contain the considered
/// lorentz vectors (must be of type ROOT::Math::PtEtaPhiMVector)
///
/// \returns a dataframe containing the new column
template <class... Lorentzvectors>
ROOT::RDF::RNode CombineP4s(ROOT::RDF::RNode df, const std::string &outputflag,
                            const Lorentzvectors &...p4s) {
    std::vector<std::string> P4List;
    utility::appendParameterPackToVector(P4List, p4s...);
    const auto np4s = sizeof...(Lorentzvectors);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        utility::PassAsVec<np4s, ROOT::Math::PtEtaPhiMVector>(
            [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p4s) {
                return Sum(p4s, ROOT::Math::PtEtaPhiMVector());
            }),
        P4List);
}

ROOT::RDF::RNode scaleP4(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvector, const float &p4_miss_sf);   
/// namespace used for mutau lorentzvectors
namespace mutau {
ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &pairname,
                       const std::vector<std::string> &muon_quantities,
                       const std::vector<std::string> &tau_quantities,
                       const std::string &muon_p4_name,
                       const std::string &tau_p4_name);
} // namespace mutau

} // namespace lorentzvectors
#endif /* GUARDLVECS_H */