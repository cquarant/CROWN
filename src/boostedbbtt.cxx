#ifndef GUARD_BOOSTEDBBTT_H
#define GUARD_BOOSTEDBBTT_H

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

}

#endif