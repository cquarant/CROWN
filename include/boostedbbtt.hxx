#ifndef GUARD_BOOSTEDBBTT_H
#define GUARD_BOOSTEDBBTT_H

namespace boostedbbtt{
    ROOT::RDF::RNode dR_tau_fatjet0(
        ROOT::RDF::RNode df, 
        const std::string &dR_tau_fatjet,
        const std::string &p0_fatjet,   
        const std::string &taupt,     
        const std::string &taueta,   
        const std::string &tauphi, 
        const std::string &taumass);
}

#endif