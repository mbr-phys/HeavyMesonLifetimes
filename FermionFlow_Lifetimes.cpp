#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

#include "Inputs.hpp"
#include "Functions.hpp"

using namespace Grid;
using namespace Hadrons;

// set fermion boundary conditions to be periodic space, antiperiodic time.
std::string boundary = "1 1 1 -1";
std::string twist = "0. 0. 0. 0.";

struct XmlPar
{
    XmlInputs::RunPar      runPar;
    XmlInputs::ConfigPar   configPar;
    XmlInputs::MdwfPar     mdwfParh;
    XmlInputs::DwfPar      dwfPars;
    XmlInputs::DwfPar      dwfParl;
    XmlInputs::MesonPar    mesonPar;
    XmlInputs::FlowPar     flowPar;
    XmlInputs::StoutPar    stoutPar;
    XmlInputs::SmearPar    smearPar;
    XmlInputs::HQEPar      hqePar;
    XmlInputs::StochasticPar stochasticPar;
};

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

    read(reader, "runPar",      xmlPar.runPar);
    read(reader, "configPar",   xmlPar.configPar);
    read(reader, "hdwfPar",     xmlPar.mdwfParh);
    read(reader, "sdwfPar",     xmlPar.dwfPars);
    read(reader, "ldwfPar",     xmlPar.dwfParl);
    read(reader, "mesonPar",    xmlPar.mesonPar);
    read(reader, "flowPar",     xmlPar.flowPar);
    read(reader, "stoutPar",    xmlPar.stoutPar);
    read(reader, "smearPar",    xmlPar.smearPar);
    read(reader, "hqePar",      xmlPar.hqePar);
    if (xmlPar.runPar.eyes == 1) {
        read(reader, "stochasticPar", xmlPar.stochasticPar);
    }

    auto fail = [](const std::string &msg) {
        std::cerr << "error: " << msg << std::endl;
        return EXIT_FAILURE;
    };

    // *****************************
    // ****** input checks *********
    // *****************************

    auto &configPar = xmlPar.configPar;
    auto &flowPar   = xmlPar.flowPar;
    auto &hqePar    = xmlPar.hqePar;
    auto &runPar    = xmlPar.runPar;
    const bool light = (xmlPar.mesonPar.light == 1);
    const bool smear = (xmlPar.smearPar.smear == 1);
    const bool eyes  = (runPar.eyes == 1);
    const bool saveWalls = (runPar.saveWalls == 1);

    const std::string chaining = xmlPar.runPar.flowChaining;
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

    if (runPar.eyes > 1 || runPar.saveWalls > 1 || runPar.wallIo32 > 1) {
        return fail("runPar.eyes, runPar.saveWalls and runPar.wallIo32 must be 0 or 1");
    }
    if (saveWalls && runPar.wallStem.empty()) {
        return fail("runPar.wallStem must be set when runPar.saveWalls = 1");
    }
    if (eyes && (xmlPar.stochasticPar.nSrcs == 0 || xmlPar.stochasticPar.nsparse == 0)) {
        return fail("stochasticPar.nSrcs and stochasticPar.nsparse must be positive when runPar.eyes = 1");
    }
    if (eyes && xmlPar.stochasticPar.eyeGammaPairs.empty()) {
        return fail("stochasticPar.eyeGammaPairs must be set (diag, all or a list of pairs) when runPar.eyes = 1");
    }

    // probe gammas (source/sink of every 3pt)
    std::vector<std::string> probes;
    std::string err = parse_probes(hqePar.probeGammas, probes);
    if (!err.empty()) {
        return fail(err);
    }
    // B* gammas (sink = source of the vector 3pts)
    std::vector<std::string> vectors;
    err = parse_vectors(hqePar.vectorGammas, vectors);
    if (!err.empty()) {
        return fail(err);
    }
    const std::string vertices   = hqePar.vertices.empty() ? "all ops" : hqePar.vertices;
    const std::string gammaPairs = hqePar.gammaPairs.empty() ? "diag" : hqePar.gammaPairs;

    // 2pt gamma pairs; an empty gamma_snk means all
    std::string gamma_pairs = "";
    if (xmlPar.mesonPar.gamma_snk == "") {
        gamma_pairs = "all";
    } else {
        std::vector<std::string> gamma_src_vec = split_gammas(xmlPar.mesonPar.gamma_src,' ');
        std::vector<std::string> gamma_snk_vec = split_gammas(xmlPar.mesonPar.gamma_snk,' ');
        // program expects same number of snk and src gammas
        if (gamma_src_vec.size() != gamma_snk_vec.size()) {
            std::cerr << "gamma_src and gamma_snk must have same number of entries:" << std::endl
                      << "    gamma_src has " << gamma_src_vec.size() << " entries," << std::endl
                      << "    gamma_snk has " << gamma_snk_vec.size() << " entries." << std::endl;

            return EXIT_FAILURE;
        }
        for (int i = 0; i < gamma_src_vec.size(); i++) {
            std::string pair = "(" + gamma_snk_vec[i] + " " + gamma_src_vec[i] + ")";
            gamma_pairs += pair;
        }
        // the probe 2pts are needed for the 3pt ratios: flowed P and A0 at the sink (PCAC, D_k>=1) for every 3pt source
        std::vector<std::string> required = {};
        for (std::string gSnk : {"Gamma5", "GammaTGamma5"}) {
            for (auto &gSrc : probes) {
                required.push_back("(" + gSnk + " " + gSrc + ")");
            }
        }
        // B* 2pts
        for (auto &g : vectors) {
            required.push_back("(" + g + " " + g + ")");
        }
        for (auto &pair : required) {
            if (gamma_pairs.find(pair) == std::string::npos) {
                return fail("mesonPar gamma pairs must contain the pair " + pair);
            }
        }
    }

    // wall pairs t_i -> t_f = (t_i + dT) mod T, looked up among the source walls
    std::vector<WallPair> pairs;
    err = make_wall_pairs(configPar, pairs);
    if (!err.empty()) {
        return fail(err);
    }

    // two-quark triads "sink vertex source"; B* triads are diagonal in the polarisation
    std::vector<std::string> triads = {};
    for (auto &v : split_gammas(vertices, ' ')) {
        for (auto &gSnk : probes) {
            for (auto &gSrc : probes) {
                triads.push_back(gSnk + " " + v + " " + gSrc);
            }
        }
        for (auto &g : vectors) {
            triads.push_back(g + " " + v + " " + g);
        }
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

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start        = configPar.begin;
    globalPar.trajCounter.end          = configPar.end;
    globalPar.trajCounter.step         = configPar.step;
    globalPar.runId                    = xmlPar.runPar.runId;
    globalPar.database.applicationDb   = xmlPar.runPar.appDb;
    globalPar.database.restoreSchedule = xmlPar.runPar.restore;
    if (xmlPar.runPar.scheduler == "naive") {
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

    // stout smear (heavy action only)
    MGauge::StoutSmearing::Par stPar;
    stPar.gauge = "trans_gauge";
    stPar.steps = xmlPar.stoutPar.steps;
    stPar.orthogDim = xmlPar.stoutPar.orthog;
    stPar.rho = xmlPar.stoutPar.rho;
    application.createModule<MGauge::StoutSmearing>("Vgauge",stPar);

    // single prec cast
    MUtilities::GaugeSinglePrecisionCast::Par fPar;
    fPar.field = "trans_gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("trans_gaugeF",fPar);

    // *****************************
    // ****** walls and props ******
    // *****************************

    std::vector<Wall> walls;
    std::vector<std::string> step0 = {};
    const std::string tag0 = std::to_string(configPar.sources[0]); // actions are named after the first wall
    for (int k = 0; k < (int)configPar.sources.size(); k++) {
        int j = (k == 0) ? 0 : 1;

        // Z2 wall, smeared for the light/strange quarks if smear == 1
        Wall w = make_wall(application, configPar.sources[k], smear, xmlPar.smearPar);

        // heavy quark
        w.q["b"] = make_moebiusdwfpropagator(application, xmlPar.mdwfParh, w.z2, j, boundary, twist, "Vgauge", "_"+w.tag);
        step0.push_back(make_wardid(application, w.q["b"], "MDWF_action_"+xmlPar.mdwfParh.quark+"_"+tag0, xmlPar.mdwfParh.mass, w.z2, "", "WardId_"+xmlPar.mdwfParh.quark+w.tag));

        // strange quark
        w.q["s"] = make_mixeddwfpropagator(application, xmlPar.dwfPars, w.gsrc, j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+w.tag);
        step0.push_back(make_wardid(application, w.q["s"], "DWF_action_"+xmlPar.dwfPars.quark+"_"+tag0, xmlPar.dwfPars.mass, w.gsrc, "", "WardId_"+xmlPar.dwfPars.quark+w.tag));

        // light quark
        if (light) {
            w.q["l"] = make_mixeddwfpropagator(application, xmlPar.dwfParl, w.gsrc, j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+w.tag);
            step0.push_back(make_wardid(application, w.q["l"], "DWF_action_"+xmlPar.dwfParl.quark+"_"+tag0, xmlPar.dwfParl.mass, w.gsrc, "", "WardId_"+xmlPar.dwfParl.quark+w.tag));
        }
        walls.push_back(w);
    }

    // flavours that are flowed, and spectators of the 3pts
    std::vector<std::string> flavs = {"b", "s"};
    std::vector<std::string> spectators = {"s"};
    if (light) {
        flavs.push_back("l");
        spectators.push_back("l");
    }

    // save the unflowed wall propagators for FermionFlow_Eyes, and the light/strange wall sources that close its
    // spectators (so FermionFlow_Eyes can use any runId for new loop noise)
    if (saveWalls) {
        if (Grid::CartesianCommunicator::RankWorld() == 0 && Hadrons::mkdir(runPar.wallStem)) {
            HADRONS_ERROR(Io, "cannot create runPar.wallStem directory '" + runPar.wallStem + "'");
        }
        for (auto &w : walls) {
            for (auto &f : flavs) {
                make_save_wall(application, w.q[f], runPar.wallStem, runPar.wallIo32 == 1);
            }
            make_save_field(application, w.gsrc, wall_source(w.t), runPar.wallStem, runPar.wallIo32 == 1);
        }
    }

    // stochastic loops for the eyes and condensates (solvers already created by the walls)
    EyeNoise noise;
    if (eyes) {
        std::map<std::string, XmlInputs::DwfPar *> loopPars = {{"s", &xmlPar.dwfPars}};
        if (light) {
            loopPars["l"] = &xmlPar.dwfParl;
        }
        noise = make_eye_noise(application, xmlPar.stochasticPar, loopPars, false, boundary, twist, "trans_gauge", "trans_gaugeF");
        make_eye_spectators(application, walls, pairs, spectators);
    }

    std::string snk_mom = "0 0 0"; // only ever zero momentum
    std::string snk_str = make_PtSnk(application, snk_mom, "");
    const std::string SM = smear ? "SP" : "PP";

    // all contractions at one flow time; label is the flow-time label, gauge the (flowed) gauge field
    auto contract = [&](const std::string &label, const std::string &gauge_tau, double tau, bool tau0, std::vector<std::string> &group) {
        // propagator of flavour f from wall w at this flow time
        auto qt = [&](const Wall &w, const std::string &f) {
            return tau0 ? w.q.at(f) : w.q.at(f) + "_t" + label;
        };

        for (auto &w : walls) {
            const std::string hq = qt(w, "b"), sq = qt(w, "s");
            const std::string suf = "_t" + w.tag + "_" + label;

            // 2pt contractions
            group.push_back(make_contraction(application, "sh_"+SM+suf, sq, hq, gamma_pairs, snk_str, "", true, false));
            group.push_back(make_contraction(application, "ss_"+SM+suf, sq, sq, gamma_pairs, snk_str, "", true, false));
            group.push_back(make_contraction(application, "hh_PP"+suf, hq, hq, gamma_pairs, snk_str, "", true, false));
            if (light) {
                const std::string lq = qt(w, "l");
                group.push_back(make_contraction(application, "lh_"+SM+suf, lq, hq, gamma_pairs, snk_str, "", true, false));
                group.push_back(make_contraction(application, "ls_"+SM+suf, lq, sq, gamma_pairs, snk_str, "", true, false));
                group.push_back(make_contraction(application, "ll_"+SM+suf, lq, lq, gamma_pairs, snk_str, "", true, false));
            }

            // 2pt contractions with smearing at source and sink, only at zero flow time
            if (tau0 && xmlPar.mesonPar.ss == 1 && smear) {
                std::string sqsm = make_SmrSrc(application, xmlPar.smearPar, sq, "trans_gauge", "_smr");
                group.push_back(make_contraction(application, "sh_SS"+suf, sqsm, hq, gamma_pairs, snk_str, "", true, false));
                group.push_back(make_contraction(application, "ss_SS"+suf, sqsm, sqsm, gamma_pairs, snk_str, "", true, false));
                if (light) {
                    std::string lqsm = make_SmrSrc(application, xmlPar.smearPar, qt(w, "l"), "trans_gauge", "_smr");
                    group.push_back(make_contraction(application, "lh_SS"+suf, lqsm, hq, gamma_pairs, snk_str, "", true, false));
                    group.push_back(make_contraction(application, "ls_SS"+suf, lqsm, sqsm, gamma_pairs, snk_str, "", true, false));
                    group.push_back(make_contraction(application, "ll_SS"+suf, lqsm, lqsm, gamma_pairs, snk_str, "", true, false));
                }
            }
        }

        // 3pt contractions for every wall pair and spectator
        for (auto &p : pairs) {
            const Wall &wi = walls[p.i], &wf = walls[p.f];
            const std::string hi = qt(wi, "b"), hf = qt(wf, "b");
            const std::string suf = "_t" + wi.tag + "_dt" + std::to_string(p.dT) + "_" + label;

            for (auto &f : spectators) {
                const std::string qi = qt(wi, f), qf = qt(wf, f);
                const std::string pre = f + "h_";

                // dim-6 weak non-eye contractions, for every probe pair
                for (auto &gIn : probes) {
                    for (auto &gOut : probes) {
                        Gamma::Algebra gI = strToVec<Gamma::Algebra>(gIn)[0], gO = strToVec<Gamma::Algebra>(gOut)[0];
                        const std::string g = "_" + probeTag(gIn) + "_" + probeTag(gOut);
                        group.push_back(make_weaknoneye(application, hi, qi, qf, hf, gI, gO, "", pre+"DeltaHs2"+g+suf));
                        group.push_back(make_weaknoneye(application, qi, hi, qf, hf, gI, gO, "", pre+"DeltaHs0_T1"+g+suf));
                        group.push_back(make_weaknoneye(application, qi, qi, hf, hf, gI, gO, "", pre+"DeltaHs0_T2"+g+suf));
                    }
                }

                // dim-6 weak non-eye contractions in the B*, diagonal in the polarisation
                for (auto &gV : vectors) {
                    Gamma::Algebra gG = strToVec<Gamma::Algebra>(gV)[0];
                    const std::string g = "_" + vectorTag(gV) + "_" + vectorTag(gV);
                    group.push_back(make_weaknoneye(application, hi, qi, qf, hf, gG, gG, "", pre+"DeltaHs2"+g+suf));
                    group.push_back(make_weaknoneye(application, qi, hi, qf, hf, gG, gG, "", pre+"DeltaHs0_T1"+g+suf));
                    group.push_back(make_weaknoneye(application, qi, qi, hf, hf, gG, gG, "", pre+"DeltaHs0_T2"+g+suf));
                }

                // two-quark HQE 3pts: flowed heavy legs, unflowed spectator, sink noise of wall t_f
                if (hqePar.twoQuark && (!tau0 || hqePar.atTau0)) {
                    group.push_back(make_hm3pt(application, pre+"HQE"+suf, gauge_tau, flowPar.bc, tau, hi, hf, wi.q.at(f), wf.gsrc, wf.t, triads));
                }

                // dim-7 DeltaB=2: DeltaHs2 legs with the Lenz-Nierste derivative specifications, one module per probe pair
                if (hqePar.dim7 && (!tau0 || hqePar.atTau0)) {
                    for (auto &gIn : probes) {
                        for (auto &gOut : probes) {
                            Gamma::Algebra gI = strToVec<Gamma::Algebra>(gIn)[0], gO = strToVec<Gamma::Algebra>(gOut)[0];
                            const std::string g = "_" + probeTag(gIn) + "_" + probeTag(gOut);
                            group.push_back(make_weaknoneye_deriv(application, hi, qi, qf, hf, gI, gO, gammaPairs, dim7Specs(), gauge_tau, flowPar.bc, pre+"DeltaHs2_dim7"+g+suf));
                        }
                    }
                }
            }
        }

        // dim-6 DeltaB=0 eyes and condensates
        if (eyes) {
            make_eyes(application, noise, walls, pairs, spectators, probes, label, gauge_tau, tau0, xmlPar.stochasticPar.c_fl, xmlPar.stochasticPar.eyeGammaPairs, group);
        }
    };

    //////////////////////////////////////
    //  Contract at zero flow time ///////
    //////////////////////////////////////
    const std::string label0 = tauLabel(0, flowPar.step_size);
    step0.push_back(make_wilsonflow0(application, "WilsonFlow_t"+label0, "trans_gauge", flowPar));
    contract(label0, "trans_gauge", 0., true, step0);

    // make a single group for all contractions at one flow time
    std::vector<std::string> corrgroups = {make_corrgroup(application, "FlowTime_t"+label0, step0)};

    /////////////////////////////////////////
    // Contract for each flow time step /////
    /////////////////////////////////////////
    // walls (and, with eyes, the noise and solutions) are flowed together
    std::vector<std::string> flowFields = {};
    for (auto &w : walls) {
        for (auto &f : flavs) {
            flowFields.push_back(w.q[f]);
        }
    }
    if (eyes) {
        std::vector<std::string> eyeFields = eye_flow_fields(noise);
        flowFields.insert(flowFields.end(), eyeFields.begin(), eyeFields.end());
    }
    std::vector<std::string> flowGroups = make_flow_chain(application, flowPar, chaining, "trans_gauge", flowFields,
        [&](const std::string &label, const std::string &gauge_tau, double tau, std::vector<std::string> &group) {
            contract(label, gauge_tau, tau, false, group);
        });
    corrgroups.insert(corrgroups.end(), flowGroups.begin(), flowGroups.end());

    // pass all groups of contractions to a single output file
    MIO::WriteResultGroup::Par WCPar;
    WCPar.results = corrgroups;
    WCPar.output = xmlPar.mesonPar.file_name;
    application.createModule<MIO::WriteResultGroup>("WriteToFile",WCPar);

    // execution
    application.saveParameterFile("RunFermionFlow_Lifetimes.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
