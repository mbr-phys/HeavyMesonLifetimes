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
    globalPar.trajCounter.start        = xmlPar.configPar.begin;
    globalPar.trajCounter.end          = xmlPar.configPar.end;
    globalPar.trajCounter.step         = xmlPar.configPar.step;
    globalPar.runId                    = xmlPar.runPar.runId;
    globalPar.database.applicationDb   = xmlPar.runPar.appDb;
    globalPar.database.restoreSchedule = xmlPar.runPar.restore;
    if (xmlPar.runPar.scheduler == "naive") {
        globalPar.scheduler.schedulerType = "naive";}
    application.setPar(globalPar);
    
    // gauge field
    MIO::LoadNersc::Par gauge;
    gauge.file = xmlPar.configPar.fileStem;
    application.createModule<MIO::LoadNersc>("gauge",gauge);

    // translat gauge field
    MGauge::Translat::Par tPar;
    tPar.gauge = "gauge";
    tPar.xvec = coordVec(xmlPar.configPar.source_loc);
    application.createModule<MGauge::Translat>("trans_gauge",tPar);

    // stout smear
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
    }

    // Z2 sources
    std::vector<std::string> sources = {};
    std::vector<std::string> srcs = {};
    for (int time : xmlPar.configPar.sources) {
        std::stringstream ss; ss << time;
        xmlPar.configPar.source_loc = "0 0 0 " + ss.str();
        sources.push_back(ss.str());
        srcs.push_back(make_Z2src(application, time, ""));
    }

    std::vector<std::string> hqs = {};
    std::vector<std::string> sqs = {};
    std::vector<std::string> lqs = {};
    std::vector<std::string> step0 = {};
    for (int k = 0; k < srcs.size(); k++) {
        int j = 0;
        if (k > 0) j = 1;
        // heavy quark
        std::string name_hq1 = make_moebiusdwfpropagator(application, xmlPar.mdwfParh, srcs[k], j, boundary, twist, "Vgauge", "_"+sources[k]);
        hqs.push_back(name_hq1);

        MContraction::WardIdentity::Par WPar;
        WPar.prop = name_hq1+"_5d";
        WPar.action = "MDWF_action_" + xmlPar.mdwfParh.quark + "_"+sources[0];
        WPar.mass = xmlPar.mdwfParh.mass;
        WPar.source = srcs[k];
        application.createModule<MContraction::WardIdentity>("WardId_"+xmlPar.mdwfParh.quark+sources[k],WPar);
        step0.push_back("WardId_"+xmlPar.mdwfParh.quark+sources[k]);

        // strange quark
        std::string name_sq1;
        std ::string gsrc;
        if (xmlPar.smearPar.smear == 1) {
            gsrc = make_SmrSrc(application, xmlPar.smearPar, srcs[k], "trans_gauge", "_"+sources[k]);
            name_sq1 = make_mixeddwfpropagator(application, xmlPar.dwfPars, gsrc, j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+sources[k]);
        } else {
            name_sq1 = make_mixeddwfpropagator(application, xmlPar.dwfPars, srcs[k], j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+sources[k]);
        }
        sqs.push_back(name_sq1);

        MContraction::WardIdentity::Par SPar;
        SPar.prop = name_sq1+"_5d";
        SPar.action = "DWF_action_" + xmlPar.dwfPars.quark + "_" + sources[0];
        SPar.mass = xmlPar.dwfPars.mass;
        if (xmlPar.smearPar.smear == 1) {
            SPar.source = gsrc;
        } else {
            SPar.source = srcs[k];
        }
        application.createModule<MContraction::WardIdentity>("WardId_"+xmlPar.dwfPars.quark+sources[k],SPar);
        step0.push_back("WardId_"+xmlPar.dwfPars.quark+sources[k]);

        if (xmlPar.mesonPar.light == 1) {
            // light quark
            std::string name_lq1;
            if (xmlPar.smearPar.smear == 1) {
                name_lq1 = make_mixeddwfpropagator(application, xmlPar.dwfParl, gsrc, j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+sources[k]);
            } else {
                name_lq1 = make_mixeddwfpropagator(application, xmlPar.dwfParl, srcs[k], j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+sources[k]);
            }
            lqs.push_back(name_lq1);

            MContraction::WardIdentity::Par LPar;
            LPar.prop = name_lq1+"_5d";
            LPar.action = "DWF_action_" + xmlPar.dwfParl.quark + "_" + sources[0];
            LPar.mass = xmlPar.dwfParl.mass;
            if (xmlPar.smearPar.smear == 1) {
                LPar.source = gsrc;
            } else {
                LPar.source = srcs[k];
            application.createModule<MContraction::WardIdentity>("WardId_"+xmlPar.dwfParl.quark+sources[k],LPar);
            step0.push_back("WardId_"+xmlPar.dwfParl.quark+sources[k]);
            }
        }
    }
    std::vector<std::string> qs = hqs; 
    qs.insert(qs.end(), sqs.begin(), sqs.end());
    if (xmlPar.mesonPar.light == 1) {
        qs.insert(qs.end(), lqs.begin(), lqs.end());
    }

    // gradient flow
    MGradientFlow::WilsonFlow::Par gfPar;
    gfPar.gauge = "trans_gauge";
    gfPar.steps = 0; 
    gfPar.step_size = xmlPar.flowPar.step_size;
    gfPar.meas_interval = 1; 

    MGradientFlow::WilsonFermionFlow::Par wfPar;
    wfPar.gauge = "trans_gauge";
//    wfPar.output = xmlPar.flowPar.output; // keep empty unless you want a separate file
    wfPar.steps = 1; 
    wfPar.step_size = xmlPar.flowPar.step_size;
    wfPar.meas_interval = xmlPar.flowPar.meas_interval;
    wfPar.props = qs;
    wfPar.bc = xmlPar.flowPar.bc;
//    wfPar.maxTau = xmlPar.flowPar.maxTau; // only needed for adaptive flow -> not implemented for fermion flow

    std::string snk_mom = "0 0 0"; // only ever zero momentum

    std::string snk_str = make_PtSnk(application, snk_mom, "");
    int stp = xmlPar.configPar.sources[1] - xmlPar.configPar.sources[0];
    int Tmax = xmlPar.configPar.Tmax;

    //////////////////////////////////////
    //  Contract at zero flow time ///////
    //////////////////////////////////////
    std::string flowMod = "FlowTime_t";
    std::string WFF = "WilsonFlow_t";
    std::string fM1(flowMod), WF1(WFF);
    // loop over quark sources
    for (int q = 0; q < hqs.size(); q++) {
        std::string name_sq1 = sqs[q], name_hq1 = hqs[q], tm = sources[q];
        
        std::string SM = "PP";
        if (xmlPar.smearPar.smear == 1) SM = "SP";
        // 2pt contractions
        step0.push_back(make_contraction(application, "sh_"+SM+"_t"+tm+"_0.00", name_sq1, name_hq1, gamma_pairs, snk_str, "", true, false));
        step0.push_back(make_contraction(application, "ss_"+SM+"_t"+tm+"_0.00", name_sq1, name_sq1, gamma_pairs, snk_str, "", true, false));
        step0.push_back(make_contraction(application, "hh_PP_t"+tm+"_0.00", name_hq1, name_hq1, gamma_pairs, snk_str, "", true, false));
        if (xmlPar.mesonPar.light == 1) {
            std::string name_lq1 = lqs[q];
            step0.push_back(make_contraction(application, "lh_"+SM+"_t"+tm+"_0.00", name_lq1, name_hq1, gamma_pairs, snk_str, "", true, false));
            step0.push_back(make_contraction(application, "ls_"+SM+"_t"+tm+"_0.00", name_lq1, name_sq1, gamma_pairs, snk_str, "", true, false));
            step0.push_back(make_contraction(application, "ll_"+SM+"_t"+tm+"_0.00", name_lq1, name_lq1, gamma_pairs, snk_str, "", true, false));
        }

        // 2pt contractions with smearing at source and sink
        // only at zero flow time
        if (xmlPar.mesonPar.ss == 1 && xmlPar.smearPar.smear == 1) {
            std::string name_sqsm = make_SmrSrc(application, xmlPar.smearPar, name_sq1, "trans_gauge", "_smr");
            step0.push_back(make_contraction(application, "sh_SS_t"+tm+"_0.00", name_sqsm, name_hq1, gamma_pairs, snk_str, "", true, false));
            step0.push_back(make_contraction(application, "ss_SS_t"+tm+"_0.00", name_sqsm, name_sqsm, gamma_pairs, snk_str, "", true, false));
            if (xmlPar.mesonPar.light == 1) {
                std::string name_lq1 = lqs[q];
                std::string name_lqsm = make_SmrSrc(application, xmlPar.smearPar, name_lq1, "trans_gauge", "_smr");
                step0.push_back(make_contraction(application, "lh_SS_t"+tm+"_0.00", name_lqsm, name_hq1, gamma_pairs, snk_str, "", true, false));
                step0.push_back(make_contraction(application, "ls_SS_t"+tm+"_0.00", name_lqsm, name_sqsm, gamma_pairs, snk_str, "", true, false));
                step0.push_back(make_contraction(application, "ll_SS_t"+tm+"_0.00", name_lqsm, name_lqsm, gamma_pairs, snk_str, "", true, false));
            }
        }

        // loop over all source-source separations for 3pt contractions
        for (int dT : xmlPar.configPar.deltaTs) {
            std::stringstream dts; dts << dT;
            int r = (stp*q + dT < Tmax) ? (q + dT/stp) : (stp*q + dT - Tmax)/stp;
            std::string name_sq2 = sqs[r], name_hq2 = hqs[r];

            // 3pt weak non-eye contractions
            step0.push_back(make_weaknoneye(application, name_hq1, name_sq1, name_sq2, name_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "sh_DeltaHs2_G5_G5_t"+tm+"_dt"+dts.str()+"_0.00"));
            step0.push_back(make_weaknoneye(application, name_sq1, name_hq1, name_sq2, name_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "sh_DeltaHs0_T1_G5_G5_t"+tm+"_dt"+dts.str()+"_0.00"));
            step0.push_back(make_weaknoneye(application, name_sq1, name_sq1, name_hq2, name_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "sh_DeltaHs0_T2_G5_G5_t"+tm+"_dt"+dts.str()+"_0.00"));
            if (xmlPar.mesonPar.light == 1) {
                std::string name_lq1 = lqs[q], name_lq2 = lqs[r];
                step0.push_back(make_weaknoneye(application, name_hq1, name_lq1, name_lq2, name_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "lh_DeltaHs2_G5_G5_t"+tm+"_dt"+dts.str()+"_0.00"));
                step0.push_back(make_weaknoneye(application, name_lq1, name_hq1, name_lq2, name_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "lh_DeltaHs0_T1_G5_G5_t"+tm+"_dt"+dts.str()+"_0.00"));
                step0.push_back(make_weaknoneye(application, name_lq1, name_lq1, name_hq2, name_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "lh_DeltaHs0_T2_G5_G5_t"+tm+"_dt"+dts.str()+"_0.00"));
            }
        }
    }
    // add gauge observables important along gradient flow
    application.createModule<MGradientFlow::WilsonFlow>(WFF+"0.00",gfPar);
    step0.push_back(WFF+"0.00");

    // make a single group for all contractions at one flow time
    std::vector<std::string> corrgroups = {make_corrgroup(application, flowMod+"0.00",step0)};

    /////////////////////////////////////////
    // Contract for each flow time step /////
    /////////////////////////////////////////
    for (int t = 1; t <= xmlPar.flowPar.steps; t++) {
        double ti = xmlPar.flowPar.step_size * t;
        //
        // after Tcoarsen, extend the measurement interval in the flow time
        if (ti == xmlPar.flowPar.Tcoarsen) xmlPar.flowPar.meas_interval *= 4;
        std::vector<std::string> qis = {};

        // pass evolved propagators to the next step in the gradient flow
        if (t != 1) {
            wfPar.gauge = WFF + "_U";
            for (int q = 0; q < qs.size(); q++) {
                std::stringstream ps; ps << q;
                qis.push_back(WFF+"_q"+ps.str()+"_1");
            }
            wfPar.props = qis;
        }
        std::stringstream st; st << std::fixed << std::setprecision(2) << ti;
        fM1 = "FlowTime_t" + st.str(); WF1 = "WilsonFlow_t" + st.str();
        std::vector<std::string> stepi = {};

        std::string time = "_t" + st.str();

        // only perform contractions for set intervals along the evolution
        // always perform contractions at final flow time
        if ((t % xmlPar.flowPar.meas_interval == 0) || (t == xmlPar.flowPar.steps)) {
            // add gauge observables important along gradient flow
            application.createModule<MGradientFlow::WilsonFermionFlow>(WF1,wfPar);
            
            // loop over quark sources
            int split = (xmlPar.mesonPar.light == 1) ? (3) : (2);
            for (int y = 0; y < qs.size()/split; y++) {
                int z = y + qs.size()/split;
                std::string flow_hq1 = qs[y]+time;
                std::string flow_sq1 = qs[z]+time;
                std::string tm = sources[y];

                // 2pt contractions
                stepi.push_back(make_contraction(application, "sh_t"+tm+"_"+st.str(), flow_sq1, flow_hq1, gamma_pairs, snk_str, "", true, false));
                stepi.push_back(make_contraction(application, "ss_t"+tm+"_"+st.str(), flow_sq1, flow_sq1, gamma_pairs, snk_str, "", true, false));
                stepi.push_back(make_contraction(application, "hh_t"+tm+"_"+st.str(), flow_hq1, flow_hq1, gamma_pairs, snk_str, "", true, false));
                if (xmlPar.mesonPar.light == 1) {
                    int zl = z + qs.size()/split;
                    std::string flow_lq1 = qs[zl]+time;
                    stepi.push_back(make_contraction(application, "lh_t"+tm+"_"+st.str(), flow_lq1, flow_hq1, gamma_pairs, snk_str, "", true, false));
                    stepi.push_back(make_contraction(application, "ls_t"+tm+"_"+st.str(), flow_lq1, flow_sq1, gamma_pairs, snk_str, "", true, false));
                    stepi.push_back(make_contraction(application, "ll_t"+tm+"_"+st.str(), flow_lq1, flow_lq1, gamma_pairs, snk_str, "", true, false));
                }

                // loop over all source-source separations for 3pt contractions
                for (int dT : xmlPar.configPar.deltaTs) {
                    std::stringstream dts; dts << dT;
                    int xy = (stp*y + dT < Tmax) ? (y + dT/stp) : (stp*y + dT - Tmax)/stp;
                    int xz = xy + qs.size()/split; 
                    std::string flow_hq2 = qs[xy]+time;
                    std::string flow_sq2 = qs[xz]+time;
                    
                    // 3pt weak non-eye contractions
                    stepi.push_back(make_weaknoneye(application, flow_hq1, flow_sq1, flow_sq2, flow_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "sh_DeltaHs2_G5_G5_t"+tm+"_dt"+dts.str()+"_"+st.str()));
                    stepi.push_back(make_weaknoneye(application, flow_sq1, flow_hq1, flow_sq2, flow_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "sh_DeltaHs0_T1_G5_G5_t"+tm+"_dt"+dts.str()+"_"+st.str()));
                    stepi.push_back(make_weaknoneye(application, flow_sq1, flow_sq1, flow_hq2, flow_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "sh_DeltaHs0_T2_G5_G5_t"+tm+"_dt"+dts.str()+"_"+st.str()));
                    if (xmlPar.mesonPar.light == 1) {
                        int zl = z + qs.size()/split;
                        std::string flow_lq1 = qs[zl]+time;
                        int xzl = xz + qs.size()/split;
                        std::string flow_lq2 = qs[xzl]+time;
                        stepi.push_back(make_weaknoneye(application, flow_hq1, flow_lq1, flow_lq2, flow_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "lh_DeltaHs2_G5_G5_t"+tm+"_dt"+dts.str()+"_"+st.str()));
                        stepi.push_back(make_weaknoneye(application, flow_lq1, flow_hq1, flow_lq2, flow_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "lh_DeltaHs0_T1_G5_G5_t"+tm+"_dt"+dts.str()+"_"+st.str()));
                        stepi.push_back(make_weaknoneye(application, flow_lq1, flow_lq1, flow_hq2, flow_hq2, Gamma::Algebra::Gamma5, Gamma::Algebra::Gamma5, "", "lh_DeltaHs0_T2_G5_G5_t"+tm+"_dt"+dts.str()+"_"+st.str()));
                    }
                }
            }
            stepi.push_back(WF1);
            
            // make a single group for all contractions at one flow time
            corrgroups.push_back(make_corrgroup(application, fM1, stepi));
        } else {
            // add gauge observables important along gradient flow
            application.createModule<MGradientFlow::WilsonFermionFlow>(WF1,wfPar);
            stepi.push_back(WF1);
            
            // make a single group for all contractions at one flow time
            corrgroups.push_back(make_corrgroup(application, fM1, stepi));
        }
        flowMod = fM1; WFF = WF1;
    }

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
