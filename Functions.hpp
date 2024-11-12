#pragma once

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "Inputs.hpp"

using namespace Grid;
using namespace Hadrons;

// ****************************
// ********* GENERAL **********
// ****************************

// check format of file stem name
std::string fileStem(std::string cname) {
    char last = cname.back();
    if (last == '/') return cname;
    else return cname+"/";
}

std::vector<std::string> split_gammas(std::string str, char del){
    std::string temp = "";
    std::vector<std::string> strings;
    for (int i = 0; i < (int)str.size(); i++) {
        if (str[i] != del) {
            temp += str[i];
        }
        else {
            strings.push_back(temp);
            temp = "";
        }
    }
    if (temp != "") {
        strings.push_back(temp);
    }
    return strings;
}

std::string removeSpaces(std::string word) {
    std::string newWord;
    for (int i = 0; i < word.length(); i++) {
        if (word[i] != ' ') {
            newWord += word[i];
        }
    }
    return newWord;
}

std::vector<unsigned int> coordVec(std::string coord) {
    std::vector<unsigned int> cV = {};
    std::vector<std::string> cVst = split_gammas(coord,' ');
    for (auto st : cVst) {
        cV.push_back(std::stoi(st));
    }
    return cV;
}

// ****************************
// ******* PROPAGATORS ********
// ****************************

// dwf propagator 
std::string make_dwfpropagator(Application &application, TestInputs::DwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string add){
    std::string solve_name = "CG_" + dwfPar.quark;
    if (i == 0) {
        // action
        MAction::DWF::Par actionPar;
        actionPar.gauge    = gauge; 
        actionPar.Ls       = dwfPar.Ls;
        actionPar.M5       = dwfPar.M5;
        actionPar.mass     = dwfPar.mass;
        actionPar.boundary = boundary;
        actionPar.twist    = twist;

        std::string act_name = "DWF_action_" + dwfPar.quark + add;
        application.createModule<MAction::DWF>(act_name, actionPar);

        // solver
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = act_name;
        solverPar.residual     = dwfPar.residual;
        solverPar.maxIteration = dwfPar.maxiter;

        application.createModule<MSolver::RBPrecCG>(solve_name,solverPar);
    }

    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = dwfPar.quark + add;
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

std::string make_mixeddwfpropagator(Application &application, TestInputs::DwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string Fgauge, std::string add){
    std::string solve_name = "MixedCGm_" + dwfPar.quark;
    if (i == 0) {
        MAction::DWFF::Par FactionPar;
        FactionPar.gauge    = Fgauge;
        FactionPar.Ls       = dwfPar.Ls;
        FactionPar.M5       = dwfPar.M5;
        FactionPar.mass     = dwfPar.mass;
        FactionPar.boundary = boundary;
        FactionPar.twist    = twist;

        std::string Fact_name = "DWFF_action_" + dwfPar.quark + add;
        application.createModule<MAction::DWFF>(Fact_name, FactionPar);

        MAction::DWF::Par actionPar;
        actionPar.gauge    = gauge;
        actionPar.Ls       = dwfPar.Ls;
        actionPar.M5       = dwfPar.M5;
        actionPar.mass     = dwfPar.mass;
        actionPar.boundary = boundary;
        actionPar.twist    = twist;

        std::string act_name = "DWF_action_" + dwfPar.quark + add;
        application.createModule<MAction::DWF>(act_name, actionPar);

        // solver
        MSolver::MixedPrecisionRBPrecCG::Par solverPar;
        solverPar.innerAction       = Fact_name;
        solverPar.outerAction       = act_name;
        solverPar.maxInnerIteration = dwfPar.maxiterIn;
        solverPar.maxOuterIteration = dwfPar.maxiterOut;
        solverPar.residual          = dwfPar.residual;

        application.createModule<MSolver::MixedPrecisionRBPrecCG>(solve_name, solverPar);
    }
    
    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = dwfPar.quark + add;
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

std::string make_moebiusdwfpropagator(Application &application, TestInputs::MdwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string add){
    std::string solve_name = "CGm_" + dwfPar.quark;
    if (i == 0) {
        MAction::MobiusDWF::Par actionPar;
        actionPar.gauge    = gauge;
        actionPar.Ls       = dwfPar.Ls;
        actionPar.M5       = dwfPar.M5;
        actionPar.mass     = dwfPar.mass;
        actionPar.b        = dwfPar.b;
        actionPar.c        = dwfPar.c;
        actionPar.boundary = boundary;
        actionPar.twist    = twist;

        std::string act_name = "MDWF_action_" + dwfPar.quark + add;
        application.createModule<MAction::MobiusDWF>(act_name, actionPar);

        // solver
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = act_name;
        solverPar.residual     = dwfPar.residual;
        solverPar.maxIteration = dwfPar.maxiter;

        application.createModule<MSolver::RBPrecCG>(solve_name, solverPar);
    }
    
    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = dwfPar.quark + add;
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

std::string make_mixedmoebiusdwfpropagator(Application &application, TestInputs::MdwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string Fgauge, std::string add){
    std::string solve_name = "MixedCGm_" + dwfPar.quark;
    if (i == 0) {
        MAction::MobiusDWFF::Par FactionPar;
        FactionPar.gauge    = Fgauge;
        FactionPar.Ls       = dwfPar.Ls;
        FactionPar.M5       = dwfPar.M5;
        FactionPar.mass     = dwfPar.mass;
        FactionPar.b        = dwfPar.b;
        FactionPar.c        = dwfPar.c;
        FactionPar.boundary = boundary;
        FactionPar.twist    = twist;

        std::string Fact_name = "MDWFF_action_" + dwfPar.quark + add;
        application.createModule<MAction::MobiusDWFF>(Fact_name, FactionPar);

        MAction::MobiusDWF::Par actionPar;
        actionPar.gauge    = gauge;
        actionPar.Ls       = dwfPar.Ls;
        actionPar.M5       = dwfPar.M5;
        actionPar.mass     = dwfPar.mass;
        actionPar.b        = dwfPar.b;
        actionPar.c        = dwfPar.c;
        actionPar.boundary = boundary;
        actionPar.twist    = twist;

        std::string act_name = "MDWF_action_" + dwfPar.quark + add;
        application.createModule<MAction::MobiusDWF>(act_name, actionPar);

        // solver
        MSolver::MixedPrecisionRBPrecCG::Par solverPar;
        solverPar.innerAction       = Fact_name;
        solverPar.outerAction       = act_name;
        solverPar.maxInnerIteration = dwfPar.maxiterIn;
        solverPar.maxOuterIteration = dwfPar.maxiterOut;
        solverPar.residual          = dwfPar.residual;

        application.createModule<MSolver::MixedPrecisionRBPrecCG>(solve_name, solverPar);
    }
    
    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = dwfPar.quark + add;
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

// ****************************
// ********* SOURCES **********
// ****************************

// Z2 wall source
std::string make_Z2src(Application &application, int time, std::string add){
    // start with standard Z2 source
    std::string z2_mom = "";
    MSource::Z2::Par z2Par;
    z2Par.tA = time; 
    z2Par.tB = time;
    std::stringstream st; st << time;
    z2_mom = "z2_p000_t" + st.str() + add;
    application.createModule<MSource::Z2>(z2_mom, z2Par);
    return z2_mom;
}

// Jacobi smeared source
std::string make_SmrSrc(Application &application, TestInputs::SmearPar &smearPar, std::string src, std::string gauge, std::string add){
    std::string sm_src = "sm_" + src + add;

    MSource::JacobiSmear::Par smr_src;
    smr_src.gauge = gauge;
    smr_src.width = smearPar.width;
    smr_src.iterations = smearPar.iterations;
    smr_src.orthog = smearPar.orthog;
    smr_src.source = src;
    application.createModule<MSource::JacobiSmear>(sm_src,smr_src);
    return sm_src;
}

// ****************************
// ********** SINKS ***********
// ****************************
std::string make_PtSnk(Application &application, std::string mom, std::string add){
    MSink::Point::Par sinkPar;
    sinkPar.mom = mom;
    std::string snk_str = "snk_p" + removeSpaces(mom) + add;
    application.createModule<MSink::ScalarPoint>(snk_str,sinkPar);
    return snk_str;
}

// ****************************
// ******* CONTRACTIONS *******
// ****************************

// general meson contraction 
std::string make_contraction(Application &application, std::string meson, std::string q1, std::string q2, std::string gammas, std::string sink, std::string fileSt, bool dump, bool snk){

    MContraction::Meson::Par contraction;
    contraction.q1 = q1;
    contraction.q2 = q2;
    contraction.gammas = gammas;
    contraction.sink   = sink;

    if (!dump) {
        contraction.output = fileSt + meson;}

    std::string name = fileSt + meson;
    if (snk) name += sink;
    application.createModule<MContraction::Meson>(name, contraction);
    return name;
}

//weaknoneye 4quark
std::string make_weaknoneye(Application &application, std::string qbl, std::string ql, std::string qbr, std::string qr, Gamma::Algebra gIn, Gamma::Algebra gOut, std::string output, std::string name) {

    MContraction::WeakNonEye3pt::Par weakPar;
    weakPar.qLeft = ql;
    weakPar.qBarLeft = qbl;
    weakPar.qRight = qr;
    weakPar.qBarRight = qbr;
    weakPar.gammaIn = gIn;
    weakPar.gammaOut = gOut;
    if (output.length() > 0) {
        weakPar.output = output;}

    application.createModule<MContraction::WeakNonEye3pt>(name,weakPar);
    return name;
}

//weakeye 4quark
std::string make_weakeye(Application &application, std::string qbl, std::string qbr, std::string qspec, std::string loop, unsigned int tOut, Gamma::Algebra gIn, Gamma::Algebra gOut, std::string output, std::string name) {

    MContraction::WeakEye3pt::Par weakPar;
    weakPar.qBarLeft = qbl;
    weakPar.qBarRight = qbr;
    weakPar.qSpectator = qspec;
    weakPar.loop = loop;
    weakPar.tOut = tOut;
    weakPar.gammaIn = gIn;
    weakPar.gammaOut = gOut;
    if (output.length() > 0) {
        weakPar.output = output;}

    application.createModule<MContraction::WeakEye3pt>(name,weakPar);
    return name;
}

//ward identity
std::string make_wardid(Application &application, std::string prop, std::string action, double mass, std::string source, std::string output, std::string name) {
    MContraction::WardIdentity::Par WPar;
    WPar.prop = prop+"_5d";
    WPar.action = action;
    WPar.mass = mass;
    WPar.source = source;
    if (output.length() > 0) {
        WPar.output = output;}

    application.createModule<MContraction::WardIdentity>(name,WPar);
    return name;
}

std::string make_corrgroup(Application &application, std::string name, std::vector<std::string> contractions) {
    MIO::CorrelatorGroup::Par CGPar;
    CGPar.contractions = contractions;

    application.createModule<MIO::CorrelatorGroup>(name,CGPar);
    return name;
}
