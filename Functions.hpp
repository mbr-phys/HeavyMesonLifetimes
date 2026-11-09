#pragma once

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "Inputs.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>

using namespace Grid;
using namespace Hadrons;

// ****************************
// ********* GENERAL **********
// ****************************

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

// probe gamma -> round-1 name token
std::string probeTag(const std::string &gamma) {
    if (gamma == "Gamma5") return "G5";
    if (gamma == "GammaTGamma5") return "GTG5";
    return "";
}

// probe gammas (source/sink of every 3pt); returns an error message, empty on success
std::string parse_probes(const std::string &probeGammas, std::vector<std::string> &probes) {
    probes = split_gammas(probeGammas, ' ');
    if (probes.empty()) {
        return "hqePar.probeGammas must be set (Gamma5 and/or GammaTGamma5)";
    }
    for (auto &g : probes) {
        if (probeTag(g).empty()) {
            return "hqePar.probeGammas: '" + g + "' is not Gamma5 or GammaTGamma5";
        }
        if (std::count(probes.begin(), probes.end(), g) > 1) {
            return "hqePar.probeGammas: '" + g + "' is repeated";
        }
    }
    return "";
}

// B* (vector) gamma -> name token
std::string vectorTag(const std::string &gamma) {
    if (gamma == "GammaX") return "GX";
    if (gamma == "GammaY") return "GY";
    if (gamma == "GammaZ") return "GZ";
    return "";
}

// B* gammas (sink = source of the vector 3pts), "none" for no B* 3pts; returns an error message, empty on success
std::string parse_vectors(const std::string &vectorGammas, std::vector<std::string> &vectors) {
    vectors.clear();
    if (vectorGammas == "none") {
        return "";
    }
    vectors = split_gammas(vectorGammas, ' ');
    if (vectors.empty()) {
        return "hqePar.vectorGammas must be set ('none' or a subset of GammaX GammaY GammaZ)";
    }
    for (auto &g : vectors) {
        if (vectorTag(g).empty()) {
            return "hqePar.vectorGammas: '" + g + "' is not GammaX, GammaY or GammaZ";
        }
        if (std::count(vectors.begin(), vectors.end(), g) > 1) {
            return "hqePar.vectorGammas: '" + g + "' is repeated";
        }
    }
    return "";
}

// ****************************
// ******* FLOW TIMES *********
// ****************************

// number of decimals needed for step_size*n to be unique for every step (at least 2, as in round 1)
int tauDigits(double step_size) {
    int d = 2;
    while (d < 8) {
        double scaled = step_size * std::pow(10., d);
        if (std::abs(scaled - std::round(scaled)) < 1e-9) break;
        d++;
    }
    return d;
}

// flow-time label of step n, used in field, module and group names
std::string tauLabel(int n, double step_size) {
    std::stringstream st; st << std::fixed << std::setprecision(tauDigits(step_size)) << step_size*n;
    return st.str();
}

// measured flow steps: every meas_interval steps, 4*meas_interval from Tcoarsen on (counted from the
// last measured step), and always the last step (as in SParameterGF)
std::vector<int> measuredSteps(const XmlInputs::FlowPar &flowPar) {
    std::vector<int> measured;
    const int base_interval   = flowPar.meas_interval;
    const int coarse_interval = 4*base_interval;
    int coarsen_step = static_cast<int>(std::ceil(flowPar.Tcoarsen / flowPar.step_size - 1e-12));
    if (coarsen_step < 1) coarsen_step = 1;
    int last = 0;
    for (int t = 1; t <= flowPar.steps; t++) {
        int interval = (t >= coarsen_step) ? coarse_interval : base_interval;
        if ((t - last >= interval) || (t == flowPar.steps)) {
            measured.push_back(t);
            last = t;
        }
    }
    return measured;
}

// ****************************
// ******* PROPAGATORS ********
// ****************************

// dwf propagator 
std::string make_dwfpropagator(Application &application, XmlInputs::DwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string add){
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

std::string make_mixeddwfpropagator(Application &application, XmlInputs::DwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string Fgauge, std::string add){
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

std::string make_moebiusdwfpropagator(Application &application, XmlInputs::MdwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string add){
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

std::string make_mixedmoebiusdwfpropagator(Application &application, XmlInputs::MdwfPar &dwfPar, std::string src, int i, std::string boundary, std::string twist, std::string gauge, std::string Fgauge, std::string add){
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
// ********** FLOW ************
// ****************************

// fermion flow of props over nSteps steps from gauge; only the last step is output, as outProps
// (also outputs <name>_U, the flowed gauge field at the last step)
std::string make_fermionflow(Application &application, std::string name, std::string gauge, std::vector<std::string> props, std::vector<std::string> outProps, int nSteps, XmlInputs::FlowPar &flowPar) {
    MGradientFlow::FermionFlow::Par wfPar;
    wfPar.gauge         = gauge;
    wfPar.props         = props;
    wfPar.defaultType   = "PropagatorField";
    wfPar.outProps      = outProps;
    wfPar.bc            = flowPar.bc;
    wfPar.steps         = nSteps;
    wfPar.step_size     = flowPar.step_size;
    wfPar.meas_interval = nSteps;
    application.createModule<MGradientFlow::FermionFlow>(name, wfPar);
    return name;
}

// gauge observables of the unflowed gauge field
std::string make_wilsonflow0(Application &application, std::string name, std::string gauge, XmlInputs::FlowPar &flowPar) {
    MGradientFlow::WilsonFlow::Par gfPar;
    gfPar.gauge          = gauge;
    gfPar.steps          = 0;
    gfPar.step_size      = flowPar.step_size;
    gfPar.meas_interval  = 1;
    gfPar.save_history   = false;
    gfPar.save_rk_stages = false;
    application.createModule<MGradientFlow::WilsonFlow>(name, gfPar);
    return name;
}

std::string make_corrgroup(Application &application, std::string name, std::vector<std::string> contractions);

// contractions at one measured flow time: (flow-time label, flowed gauge field, flow time, result group)
typedef std::function<void(const std::string &, const std::string &, double, std::vector<std::string> &)> FlowContraction;

// flow fields from gauge along the flow steps and contract at the measured steps; returns one result group per FermionFlow module.
// The flowed fields at step n are <field>_t<label(n)>, the flowed gauge field <WilsonFlow_t<label(n)>>_U.
// chaining = "step":     one FermionFlow module per step
//          = "interval": one FermionFlow module per measured step, covering the steps since the last one
std::vector<std::string> make_flow_chain(Application &application, XmlInputs::FlowPar &flowPar, const std::string &chaining, const std::string &gauge, const std::vector<std::string> &fields, const FlowContraction &contract) {
    std::vector<std::string> groups;
    std::vector<int> measured = measuredSteps(flowPar);
    std::vector<std::string> prevProps = fields;
    std::string prevGauge = gauge;
    int prevStep = 0;

    for (int n = 1; n <= flowPar.steps; n++) {
        const bool meas = (std::find(measured.begin(), measured.end(), n) != measured.end());
        if (chaining == "interval" && !meas) continue;

        const std::string label = tauLabel(n, flowPar.step_size);
        const std::string flowName = "WilsonFlow_t" + label;
        std::vector<std::string> outProps = {};
        for (auto &f : fields) {
            outProps.push_back(f + "_t" + label);
        }

        // flowed fields and gauge observables; <flowName>_U is the flowed gauge field at this step
        std::vector<std::string> stepi = {make_fermionflow(application, flowName, prevGauge, prevProps, outProps, n - prevStep, flowPar)};

        // only perform contractions for the measured steps (always includes the final flow time)
        if (meas) {
            contract(label, flowName + "_U", flowPar.step_size*n, stepi);
        }

        // make a single group for all contractions at one flow time
        groups.push_back(make_corrgroup(application, "FlowTime_t" + label, stepi));

        prevProps = outProps;
        prevGauge = flowName + "_U";
        prevStep  = n;
    }
    return groups;
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
std::string make_SmrSrc(Application &application, XmlInputs::SmearPar &smearPar, std::string src, std::string gauge, std::string add){
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
// ********** WALLS ***********
// ****************************

// one source wall: its time, sources and (unflowed) propagators
struct Wall
{
    int t;
    std::string tag;   // t as a string
    std::string z2;    // Z2 wall source (heavy quark source)
    std::string gsrc;  // light/strange source: smeared Z2 if smear == 1, else z2
    std::map<std::string, std::string> q; // "b", "s", "l" -> propagator name
};

// 3pt wall pair t_i -> t_f = (t_i + dT) mod T
struct WallPair
{
    int i, f, dT;
};

// wall at time t with its sources (Z2, smeared for the light/strange quarks if smear); no propagators.
// The Z2 noise is seeded from runId, module name and trajectory, so the same runId gives the same walls in every executable.
Wall make_wall(Application &application, int t, bool smear, XmlInputs::SmearPar &smearPar) {
    Wall w;
    w.t    = t;
    w.tag  = std::to_string(t);
    w.z2   = make_Z2src(application, t, "");
    w.gsrc = smear ? make_SmrSrc(application, smearPar, w.z2, "trans_gauge", "_"+w.tag) : w.z2;
    return w;
}

// wall pairs t_i -> t_f = (t_i + dT) mod T, looked up among the source walls; returns an error message, empty on success
std::string make_wall_pairs(const XmlInputs::ConfigPar &configPar, std::vector<WallPair> &pairs) {
    pairs.clear();
    for (int i = 0; i < (int)configPar.sources.size(); i++) {
        for (int dT : configPar.deltaTs) {
            int tf = (configPar.sources[i] + dT) % (int)configPar.Tmax;
            auto it = std::find(configPar.sources.begin(), configPar.sources.end(), tf);
            if (it == configPar.sources.end()) {
                return "no source wall at t_f = " + std::to_string(tf) + " (t_i = "
                       + std::to_string(configPar.sources[i]) + ", dT = " + std::to_string(dT) + ")";
            }
            pairs.push_back({i, (int)std::distance(configPar.sources.begin(), it), dT});
        }
    }
    return "";
}

// file stem of a saved wall propagator
std::string wall_file(const std::string &wallStem, const std::string &prop) {
    return wallStem + "/" + prop;
}

// file name of the saved source of wall t: the field the light/strange walls were solved on (smeared Z2 if smear == 1,
// else Z2), against which FermionFlow_Eyes closes the spectators. Saved so that the eyes do not depend on the runId seed.
std::string wall_source(int t) {
    return "src_" + std::to_string(t);
}

// save an (unflowed) 4D PropagatorField to <wallStem>/<file>.<traj>.bin
std::string make_save_field(Application &application, const std::string &field, const std::string &file, const std::string &wallStem, bool io32) {
    MIO::SavePropagator::Par savePar;
    savePar.name     = field;
    savePar.fileStem = wall_file(wallStem, file);
    std::string name = "save_" + field;
    if (io32) {
        application.createModule<MIO::SavePropagatorIo32>(name, savePar);
    } else {
        application.createModule<MIO::SavePropagator>(name, savePar);
    }
    return name;
}

// load a field saved by make_save_field; the loaded field is called field
std::string make_load_field(Application &application, const std::string &field, const std::string &file, const std::string &wallStem, bool io32) {
    MIO::LoadPropagator::Par loadPar;
    loadPar.name     = field;
    loadPar.Ls       = 1;
    loadPar.fileStem = wall_file(wallStem, file);
    std::string name = "load_" + field;
    if (io32) {
        application.createModule<MIO::LoadPropagatorIo32>(name, loadPar);
    } else {
        application.createModule<MIO::LoadPropagator>(name, loadPar);
    }
    return field;
}

// save an (unflowed) wall propagator to <wallStem>/<prop>.<traj>.bin
std::string make_save_wall(Application &application, const std::string &prop, const std::string &wallStem, bool io32) {
    return make_save_field(application, prop, prop, wallStem, io32);
}

// load a wall propagator saved by make_save_wall; the loaded field is called prop
std::string make_load_wall(Application &application, const std::string &prop, const std::string &wallStem, bool io32) {
    return make_load_field(application, prop, prop, wallStem, io32);
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

//weakeye 4quark
std::string make_weakeye(Application &application, std::string qbl, std::string qbr, std::string qspec, std::string loop, unsigned int tOut, Gamma::Algebra gIn, Gamma::Algebra gOut, std::string gammaPairs, std::string output, std::string name) {

    MContraction::WeakEye3pt::Par weakPar;
    weakPar.qBarLeft = qbl;
    weakPar.qBarRight = qbr;
    weakPar.qSpectator = qspec;
    weakPar.loop = loop;
    weakPar.tOut = tOut;
    weakPar.gammaIn = gIn;
    weakPar.gammaOut = gOut;
    weakPar.gammaPairs = gammaPairs;
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

// HeavyMeson3pt: two-quark (HQE) 3pts, zero momentum
std::string make_hm3pt(Application &application, std::string name, std::string gauge, int bc, double flowTime, std::string qSrc, std::string qSnk, std::string spectator, std::string sinkNoise, unsigned int tSnk, std::vector<std::string> triads) {
    MContraction::HeavyMeson3pt::Par hmPar;
    hmPar.gauge     = gauge;
    hmPar.bc        = bc;
    hmPar.flowTime  = flowTime;
    hmPar.qSrc      = qSrc;
    hmPar.qSnk      = qSnk;
    hmPar.spectator = spectator;
    hmPar.sinkNoise = sinkNoise;
    hmPar.tSnk      = tSnk;
    hmPar.gamma     = triads;
    application.createModule<MContraction::HeavyMeson3pt>(name, hmPar);
    return name;
}

// derivative specification of a WeakNonEye3ptDeriv contraction
MContraction::WeakDerivSpec deriv(std::string legs, std::string dirs, std::string couple) {
    MContraction::WeakDerivSpec d;
    d.legs   = legs;
    d.dirs   = dirs;
    d.couple = couple;
    return d;
}

// dim-7 DeltaB=2 derivative specifications, Beneke-Buchalla-Dunietz basis R_0..R_4, R~_1..R~_4
// (hep-ph/9605259; HQE_plan.md section 3.2)
// legs: qBarLeft = h_i, qLeft = q_i, qBarRight = q_f, qRight = h_f
// A-D are split into time (3) and space (0 1 2); with D_0 = iD_4 the Minkowski D_rho...D^rho of R_2, R_3
// is minus the sum of the two Euclidean parts.
std::vector<MContraction::WeakDerivSpec> dim7Specs(void) {
    return {
        deriv("none", "", ""),                                // dim 6: R_0, R_1, R~_1
        deriv("qLeft", "3", "none"),                          // leading-power R_2, R_3 (both Wick assignments)
        deriv("qBarRight", "3", "none"),
        deriv("qBarLeft", "3", "none"),                       // D_4 on each bbar leg: -/+ M check
        deriv("qRight", "3", "none"),
        deriv("qBarLeft qLeft", "3", "none"),                 // A: full R_2, R_3, R~_2, R~_3
        deriv("qBarLeft qLeft", "0 1 2", "none"),
        deriv("qRight qLeft", "3", "none"),                   // B
        deriv("qRight qLeft", "0 1 2", "none"),
        deriv("qBarRight qRight", "3", "none"),               // C
        deriv("qBarRight qRight", "0 1 2", "none"),
        deriv("qBarLeft qBarRight", "3", "none"),             // D
        deriv("qBarLeft qBarRight", "0 1 2", "none"),
        // R_4, R~_4: D_mu on s, gamma_mu in the other bilinear. With s = qLeft, G1 sits in the s bilinear in
        // trace 2 (R_4: couple G2) and in the other bilinear in trace 1 (R~_4: couple G1); with s = qBarRight
        // the other way round. Each spec's other trace is a Fierz partner.
        deriv("qLeft", "0 1 2 3", "G2"),
        deriv("qLeft", "0 1 2 3", "G1"),
        deriv("qBarRight", "0 1 2 3", "G1"),
        deriv("qBarRight", "0 1 2 3", "G2")
    };
}

// four-quark operator gamma pairs; empty or "hqe" = the 16 diagonal G x G (O_1, O_2, T_1, T_2 and the scalar-class
// penguins) plus the 16 parity-even (G, gamma_4 G) pairs that the Fierzed colour-octet vector-class penguins
// (bbar t^a b)(qbar gamma_4 t^a q) need, 12 of them with a sigma; anything else (diag, all, a list) is passed through
std::string gamma_pairs(const std::string &pairs) {
    if (!pairs.empty() && pairs != "hqe") {
        return pairs;
    }
    const std::vector<std::string> diag = {
        "Identity", "Gamma5", "GammaX", "GammaY", "GammaZ", "GammaT",
        "GammaXGamma5", "GammaYGamma5", "GammaZGamma5", "GammaTGamma5",
        "SigmaXY", "SigmaXZ", "SigmaYZ", "SigmaXT", "SigmaYT", "SigmaZT"
    };
    const std::vector<std::pair<std::string, std::string>> vectorClass = {
        {"Identity", "GammaT"}, {"Gamma5", "GammaTGamma5"},
        {"GammaX", "SigmaXT"}, {"GammaY", "SigmaYT"}, {"GammaZ", "SigmaZT"},
        {"GammaXGamma5", "SigmaYZ"}, {"GammaYGamma5", "SigmaXZ"}, {"GammaZGamma5", "SigmaXY"}
    };
    std::string list;
    for (auto &g : diag) {
        list += "(" + g + " " + g + ")";
    }
    for (auto &p : vectorClass) {
        list += "(" + p.first + " " + p.second + ")(" + p.second + " " + p.first + ")";
    }
    return list;
}

// WeakNonEye3ptDeriv: four-quark 3pts with gamma pairs and derivatives on the legs, zero momentum
std::string make_weaknoneye_deriv(Application &application, std::string qbl, std::string ql, std::string qbr, std::string qr, Gamma::Algebra gIn, Gamma::Algebra gOut, std::string gammaPairs, std::vector<MContraction::WeakDerivSpec> specs, std::string gauge, int bc, std::string name) {
    MContraction::WeakNonEye3ptDeriv::Par weakPar;
    weakPar.qBarLeft    = qbl;
    weakPar.qLeft       = ql;
    weakPar.qBarRight   = qbr;
    weakPar.qRight      = qr;
    weakPar.gammaIn     = gIn;
    weakPar.gammaOut    = gOut;
    weakPar.gammaPairs  = gammaPairs;
    weakPar.derivatives = specs;
    weakPar.gauge       = gauge;
    weakPar.bc          = bc;
    application.createModule<MContraction::WeakNonEye3ptDeriv>(name, weakPar);
    return name;
}

std::string make_corrgroup(Application &application, std::string name, std::vector<std::string> contractions) {
    MIO::ResultGroup::Par CGPar;
    CGPar.results = contractions;

    application.createModule<MIO::ResultGroup>(name,CGPar);
    return name;
}

// ****************************
// *********** EYES ***********
// ****************************
// Dim-6 DeltaB=0 eye diagrams (WeakEye3pt) and the SParameterGF condensates, on the same stochastic loops (HQE_plan.md section 7.1).

// stochastic noise components and their solutions
struct EyeNoise
{
    unsigned int nSrcs;                                     // noise hits
    std::vector<std::string> eta;                           // diluted noise components (nSrcs*nsparse^4)
    std::map<std::string, std::vector<std::string>> phi;    // loop flavour -> D^-1 eta, one per component
};

// sparse spin-colour diagonal Z2 noise, diluted and unpacked into nsrc*nsparse^4 PropagatorField components
// (as in SParameterGF; the estimators sum all components and divide by nsrc)
std::vector<std::string> make_stochastic_noise_diluted(Application &application, const std::string &base_name, unsigned int nsrc, unsigned int nsparse) {
    unsigned int nDilutions = nsrc * std::pow(nsparse, 4);

    MNoise::SparseSpinColorDiagonal::Par sparsePar;
    sparsePar.nsrc    = nsrc;
    sparsePar.nsparse = nsparse;
    std::string sparse_name = base_name + "_sparse";
    application.createModule<MNoise::SparseSpinColorDiagonal>(sparse_name, sparsePar);

    MSource::Z2Diluted::Par dilutedPar;
    dilutedPar.noise = sparse_name;
    application.createModule<MSource::Z2Diluted>(base_name, dilutedPar);

    MUtilities::PropagatorVectorUnpack::Par unpackPar;
    unpackPar.input = base_name;
    unpackPar.size  = nDilutions;
    std::string unpack_name = base_name + "_unpacked";
    application.createModule<MUtilities::PropagatorVectorUnpack>(unpack_name, unpackPar);

    std::vector<std::string> sources;
    for (unsigned int i = 0; i < nDilutions; i++) {
        sources.push_back(unpack_name + "_" + std::to_string(i));
    }
    return sources;
}

// noise and one solve per component and loop flavour (DWF, mixed precision, as the wall propagators).
// newSolvers: create the actions and solvers (false if the wall propagators already created them)
EyeNoise make_eye_noise(Application &application, XmlInputs::StochasticPar &stochPar, std::map<std::string, XmlInputs::DwfPar *> dwfPars, bool newSolvers, std::string boundary, std::string twist, std::string gauge, std::string gaugeF) {
    EyeNoise noise;
    noise.nSrcs = stochPar.nSrcs;
    noise.eta   = make_stochastic_noise_diluted(application, "eta", stochPar.nSrcs, stochPar.nsparse);
    for (auto &fp : dwfPars) {
        for (unsigned int k = 0; k < noise.eta.size(); k++) {
            int j = (newSolvers && k == 0) ? 0 : 1;
            noise.phi[fp.first].push_back(make_mixeddwfpropagator(application, *fp.second, noise.eta[k], j, boundary, twist, gauge, gaugeF, "_stoch_" + std::to_string(k)));
        }
    }
    return noise;
}

// fields to flow with the walls: noise and solutions (flowed together so that the loop is K D^-1 K^dagger)
std::vector<std::string> eye_flow_fields(const EyeNoise &noise) {
    std::vector<std::string> fields = noise.eta;
    for (auto &fp : noise.phi) {
        fields.insert(fields.end(), fp.second.begin(), fp.second.end());
    }
    return fields;
}

// Dslash of a field (for the <eta^dagger Dslash phi> condensate)
void make_dslash_field(Application &application, const std::string &input_field, const std::string &gauge_field, const std::string &output_name) {
    MUtilities::DslashFieldPropagator::Par dslashPar;
    dslashPar.input  = input_field;
    dslashPar.gauge  = gauge_field;
    dslashPar.output = output_name;
    application.createModule<MUtilities::DslashFieldPropagator>("Dslash_" + input_field, dslashPar);
}

// stochastic condensate -<eta^dagger Gamma phi> (+ c_fl <eta^dagger eta> for Gamma = Identity)
std::string make_stochastic_condensate(Application &application, const std::string &moduleName, const std::string &eta_field, const std::string &phi_field, const std::string &gammas, double c_fl, const std::string &output_name, bool dump) {
    MContraction::StochasticCondensatePropagator::Par condPar;
    condPar.eta    = eta_field;
    condPar.phi    = phi_field;
    condPar.gammas = gammas;
    condPar.c_fl   = c_fl;
    if (!dump) {
        condPar.output = output_name;
    }
    application.createModule<MContraction::StochasticCondensatePropagator>(moduleName, condPar);
    return moduleName;
}

// loop q adj(eta) of one noise component
std::string make_loop(Application &application, const std::string &name, const std::string &q, const std::string &eta) {
    MContraction::Loop::Par loopPar;
    loopPar.q   = q;
    loopPar.eta = eta;
    application.createModule<MContraction::Loop>(name, loopPar);
    return name;
}

// sum_i coeffs[i] fields[i]
std::string make_linearcomb(Application &application, const std::string &name, const std::vector<std::string> &fields, const std::vector<double> &coeffs) {
    MUtilities::LinearComb::Par lcPar;
    lcPar.fields = fields;
    lcPar.coeffs = coeffs;
    application.createModule<MUtilities::LinearComb>(name, lcPar);
    return name;
}

// spectator f from wall t_i closed at wall t_f against adj(sink noise of t_f) (SlicedPropagator)
std::string eye_spectator(const Wall &wi, const Wall &wf, const std::string &f) {
    return f + "_t" + wi.tag + "_snk" + wf.tag;
}

// eye spectators for every wall pair and spectator flavour, from the unflowed spectator propagators; flow-time independent
void make_eye_spectators(Application &application, const std::vector<Wall> &walls, const std::vector<WallPair> &pairs, const std::vector<std::string> &spectators) {
    auto &vm = VirtualMachine::getInstance();
    for (auto &p : pairs) {
        const Wall &wi = walls[p.i], &wf = walls[p.f];
        const std::string sink = "adjsnk_" + wf.gsrc;
        if (!vm.hasModule(sink)) {
            MSink::SourceSink::Par sinkPar;
            sinkPar.source   = wf.gsrc;
            sinkPar.sourceOp = "adj";
            application.createModule<MSink::SourceSink>(sink, sinkPar);
        }
        for (auto &f : spectators) {
            const std::string name = eye_spectator(wi, wf, f);
            if (!vm.hasModule(name)) {
                MSink::Smear::Par smearPar;
                smearPar.q    = wi.q.at(f);
                smearPar.sink = sink;
                application.createModule<MSink::Smear>(name, smearPar);
            }
        }
    }
}

// condensates, loops and eye contractions at one flow time (label; tau0: unflowed fields), all added to group.
// Fields at flow time label are <field>_t<label> (make_flow_chain); gauge is the gauge field at this flow time.
// Names (suffix _t<t_i>_dt<dT>_<label> as for the other 3pts):
//   condensate_<g>_stoch_<k>_t<label>, dslash_condensate_<g>_stoch_<k>_t<label>     (as in SParameterGF)
//   <f>h_DeltaHs0_eye_<g>_<G>_<G>_t<t_i>_dt<dT>_<label>                             (WeakEye3pt, loop flavour g, gamma pairs gammaPairs)
// <G>_<G> runs over every pair of probes (pseudoscalar walls) and over the diagonal B* walls <V>_<V> (vectors).
void make_eyes(Application &application, const EyeNoise &noise, const std::vector<Wall> &walls, const std::vector<WallPair> &pairs, const std::vector<std::string> &spectators, const std::vector<std::string> &probes, const std::vector<std::string> &vectors, const std::string &label, const std::string &gauge, bool tau0, double c_fl, const std::string &gammaPairs, std::vector<std::string> &group) {
    auto at = [&](const std::string &field) {
        return tau0 ? field : field + "_t" + label;
    };
    const std::string tl = "_t" + label;
    std::map<std::string, std::string> loop;

    for (auto &fp : noise.phi) {
        const std::string &g = fp.first;
        std::vector<std::string> loops;
        for (unsigned int k = 0; k < noise.eta.size(); k++) {
            const std::string eta = at(noise.eta[k]), phi = at(fp.second[k]);
            const std::string sk  = g + "_stoch_" + std::to_string(k) + tl;

            // condensates: all gammas, and the Dslash (Z_chi) trace; c_fl only at nonzero flow time, as in SParameterGF
            group.push_back(make_stochastic_condensate(application, "condensate_" + sk, eta, phi, "all", tau0 ? 0. : c_fl, "", true));
            make_dslash_field(application, phi, gauge, "dslash_phi_" + sk);
            group.push_back(make_stochastic_condensate(application, "dslash_condensate_" + sk, eta, "dslash_phi_" + sk, "Identity", 0., "", true));

            loops.push_back(make_loop(application, "loop_" + sk, phi, eta));
        }
        // loop estimate: sum over the components, divided by the number of hits
        loop[g] = make_linearcomb(application, "loop_" + g + tl, loops, std::vector<double>(loops.size(), 1./noise.nSrcs));
    }

    for (auto &p : pairs) {
        const Wall &wi = walls[p.i], &wf = walls[p.f];
        const std::string hi = at(wi.q.at("b")), hf = at(wf.q.at("b"));
        const std::string suf = "_t" + wi.tag + "_dt" + std::to_string(p.dT) + "_" + label;

        for (auto &f : spectators) {
            for (auto &gl : loop) {
                for (auto &gIn : probes) {
                    for (auto &gOut : probes) {
                        Gamma::Algebra gI = strToVec<Gamma::Algebra>(gIn)[0], gO = strToVec<Gamma::Algebra>(gOut)[0];
                        const std::string g = "_" + probeTag(gIn) + "_" + probeTag(gOut);
                        group.push_back(make_weakeye(application, hi, hf, eye_spectator(wi, wf, f), gl.second, wf.t, gI, gO, gammaPairs, "", f + "h_DeltaHs0_eye_" + gl.first + g + suf));
                    }
                }
                // B*, diagonal in the polarisation
                for (auto &gV : vectors) {
                    Gamma::Algebra gG = strToVec<Gamma::Algebra>(gV)[0];
                    const std::string g = "_" + vectorTag(gV) + "_" + vectorTag(gV);
                    group.push_back(make_weakeye(application, hi, hf, eye_spectator(wi, wf, f), gl.second, wf.t, gG, gG, gammaPairs, "", f + "h_DeltaHs0_eye_" + gl.first + g + suf));
                }
            }
        }
    }
}
