#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

#include "Inputs.hpp"
#include "Functions.hpp"

using namespace Grid;
using namespace Hadrons;

// Dim-6 DeltaB=0 eye diagrams and condensates only, from the unflowed wall propagators saved by
// FermionFlow_Lifetimes (runPar.saveWalls = 1). Uses the same input file: runPar.runId must be the same as
// in the run that saved the walls, so that the Z2 wall noise (needed for the spectator sink) is identical.
// The stochastic noise, solves, flow and contractions are the same as with runPar.eyes = 1 in FermionFlow_Lifetimes.

// set fermion boundary conditions to be periodic space, antiperiodic time.
std::string boundary = "1 1 1 -1";
std::string twist = "0. 0. 0. 0.";

struct XmlPar
{
    XmlInputs::RunPar        runPar;
    XmlInputs::ConfigPar     configPar;
    XmlInputs::MdwfPar       mdwfParh;
    XmlInputs::DwfPar        dwfPars;
    XmlInputs::DwfPar        dwfParl;
    XmlInputs::MesonPar      mesonPar;
    XmlInputs::FlowPar       flowPar;
    XmlInputs::SmearPar      smearPar;
    XmlInputs::HQEPar        hqePar;
    XmlInputs::StochasticPar stochasticPar;
};

// <stem>.db -> <stem>_eyes.db, else <stem>_eyes (empty stays empty)
std::string eyesDb(const std::string &db) {
    if (db.empty()) return db;
    const std::string ext = ".db";
    if (db.size() > ext.size() && db.compare(db.size() - ext.size(), ext.size(), ext) == 0) {
        return db.substr(0, db.size() - ext.size()) + "_eyes" + ext;
    }
    return db + "_eyes";
}

int main(int argc, char *argv[])
{
    //program expects parameter file as first argument of command line
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;

        return EXIT_FAILURE;
    }
    std::string parFilename;
    parFilename = argv[1];

    XmlPar xmlPar;
    XmlReader reader(parFilename);

    read(reader, "runPar",        xmlPar.runPar);
    read(reader, "configPar",     xmlPar.configPar);
    read(reader, "hdwfPar",       xmlPar.mdwfParh);
    read(reader, "sdwfPar",       xmlPar.dwfPars);
    read(reader, "ldwfPar",       xmlPar.dwfParl);
    read(reader, "mesonPar",      xmlPar.mesonPar);
    read(reader, "flowPar",       xmlPar.flowPar);
    read(reader, "smearPar",      xmlPar.smearPar);
    read(reader, "hqePar",        xmlPar.hqePar);
    read(reader, "stochasticPar", xmlPar.stochasticPar);

    auto fail = [](const std::string &msg) {
        std::cerr << "error: " << msg << std::endl;
        return EXIT_FAILURE;
    };

    // *****************************
    // ****** input checks *********
    // *****************************

    auto &configPar = xmlPar.configPar;
    auto &flowPar   = xmlPar.flowPar;
    auto &runPar    = xmlPar.runPar;
    const bool light = (xmlPar.mesonPar.light == 1);
    const bool smear = (xmlPar.smearPar.smear == 1);

    const std::string chaining = runPar.flowChaining;
    if (chaining != "step" && chaining != "interval") {
        return fail("runPar.flowChaining must be 'step' or 'interval'");
    }
    if (flowPar.steps < 0) {
        return fail("flowPar.steps must be non-negative");
    }
    if (flowPar.steps > 0 && (flowPar.step_size <= 0. || flowPar.meas_interval < 1)) {
        return fail("flowPar.step_size must be positive and flowPar.meas_interval at least 1");
    }
    if (flowPar.bc != 1 && flowPar.bc != -1) {
        return fail("flowPar.bc must be +1 or -1");
    }
    if (configPar.sources.empty() || configPar.Tmax == 0) {
        return fail("configPar.sources must not be empty and configPar.Tmax must be set");
    }
    for (int t : configPar.sources) {
        if (t < 0 || t >= (int)configPar.Tmax) {
            return fail("source time " + std::to_string(t) + " outside [0, Tmax)");
        }
    }
    if (runPar.wallStem.empty() || runPar.wallIo32 > 1) {
        return fail("runPar.wallStem must be set and runPar.wallIo32 must be 0 or 1");
    }
    if (xmlPar.stochasticPar.nSrcs == 0 || xmlPar.stochasticPar.nsparse == 0) {
        return fail("stochasticPar.nSrcs and stochasticPar.nsparse must be positive");
    }
    if (xmlPar.stochasticPar.eyeGammaPairs.empty()) {
        return fail("stochasticPar.eyeGammaPairs must be set (diag, all or a list of pairs)");
    }

    std::vector<std::string> probes;
    std::string err = parse_probes(xmlPar.hqePar.probeGammas, probes);
    if (!err.empty()) {
        return fail(err);
    }

    std::vector<WallPair> pairs;
    err = make_wall_pairs(configPar, pairs);
    if (!err.empty()) {
        return fail(err);
    }

    // *****************************
    // ****** initialization *******
    // *****************************

    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    // run setup
    Application application;

    // global parameters (same runId as the run that saved the walls)
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start        = configPar.begin;
    globalPar.trajCounter.end          = configPar.end;
    globalPar.trajCounter.step         = configPar.step;
    globalPar.runId                    = runPar.runId;
    globalPar.database.applicationDb   = eyesDb(runPar.appDb);
    globalPar.database.restoreSchedule = runPar.restore;
    if (runPar.scheduler == "naive") {
        globalPar.scheduler.schedulerType = "naive";}
    application.setPar(globalPar);

    // gauge field
    MIO::LoadNersc::Par gauge;
    gauge.file = configPar.fileStem;
    application.createModule<MIO::LoadNersc>("gauge",gauge);

    // translat gauge field
    MUtilities::CShiftGaugeField::Par tPar;
    tPar.q = "gauge";
    tPar.shift = configPar.source_loc;
    application.createModule<MUtilities::CShiftGaugeField>("trans_gauge",tPar);

    // single prec cast
    MUtilities::GaugeSinglePrecisionCast::Par fPar;
    fPar.field = "trans_gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("trans_gaugeF",fPar);

    // *****************************
    // ****** walls and props ******
    // *****************************

    // flavours loaded from disk: heavy (flowed here) and spectators (unflowed)
    std::vector<std::string> spectators = {"s"};
    if (light) {
        spectators.push_back("l");
    }
    std::vector<std::string> flavs = {"b"};
    flavs.insert(flavs.end(), spectators.begin(), spectators.end());

    // same sources as in FermionFlow_Lifetimes; propagator names as there (<quark>_<t>)
    std::map<std::string, std::string> quark = {{"b", xmlPar.mdwfParh.quark}, {"s", xmlPar.dwfPars.quark}, {"l", xmlPar.dwfParl.quark}};
    std::vector<Wall> walls;
    for (int t : configPar.sources) {
        Wall w = make_wall(application, t, smear, xmlPar.smearPar);
        for (auto &f : flavs) {
            w.q[f] = make_load_wall(application, quark[f] + "_" + w.tag, runPar.wallStem, runPar.wallIo32 == 1);
        }
        walls.push_back(w);
    }

    // stochastic loops (this executable creates the solvers)
    std::map<std::string, XmlInputs::DwfPar *> loopPars = {{"s", &xmlPar.dwfPars}};
    if (light) {
        loopPars["l"] = &xmlPar.dwfParl;
    }
    EyeNoise noise = make_eye_noise(application, xmlPar.stochasticPar, loopPars, true, boundary, twist, "trans_gauge", "trans_gaugeF");
    make_eye_spectators(application, walls, pairs, spectators);

    //////////////////////////////////////
    //  Contract at zero flow time ///////
    //////////////////////////////////////
    const std::string label0 = tauLabel(0, flowPar.step_size);
    std::vector<std::string> step0 = {make_wilsonflow0(application, "WilsonFlow_t"+label0, "trans_gauge", flowPar)};
    make_eyes(application, noise, walls, pairs, spectators, probes, label0, "trans_gauge", true, xmlPar.stochasticPar.c_fl, xmlPar.stochasticPar.eyeGammaPairs, step0);
    std::vector<std::string> corrgroups = {make_corrgroup(application, "FlowTime_t"+label0, step0)};

    /////////////////////////////////////////
    // Contract for each flow time step /////
    /////////////////////////////////////////
    // heavy walls, noise and solutions are flowed together
    std::vector<std::string> flowFields = {};
    for (auto &w : walls) {
        flowFields.push_back(w.q["b"]);
    }
    std::vector<std::string> eyeFields = eye_flow_fields(noise);
    flowFields.insert(flowFields.end(), eyeFields.begin(), eyeFields.end());
    std::vector<std::string> flowGroups = make_flow_chain(application, flowPar, chaining, "trans_gauge", flowFields,
        [&](const std::string &label, const std::string &gauge_tau, double tau, std::vector<std::string> &group) {
            make_eyes(application, noise, walls, pairs, spectators, probes, label, gauge_tau, false, xmlPar.stochasticPar.c_fl, xmlPar.stochasticPar.eyeGammaPairs, group);
        });
    corrgroups.insert(corrgroups.end(), flowGroups.begin(), flowGroups.end());

    // pass all groups of contractions to a single output file
    MIO::WriteResultGroup::Par WCPar;
    WCPar.results = corrgroups;
    WCPar.output = xmlPar.mesonPar.file_name + "_eyes";
    application.createModule<MIO::WriteResultGroup>("WriteToFile",WCPar);

    // execution
    application.saveParameterFile("RunFermionFlow_Eyes.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
