#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

#include <cmath>

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
    XmlInputs::MdwfPar     mdwfPar1;
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
    read(reader, "dwfPars",     xmlPar.mdwfPars);
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
    MUtilities::CShiftGaugeField::Par tPar;
    tPar.q = "gauge";
    tPar.shift = xmlPar.configPar.source_loc;
    application.createModule<MUtilities::CShiftGaugeField>("trans_gauge",tPar);

    // stout smear
    MGauge::StoutSmearing::Par stPar;
    stPar.gauge = "gauge";
    stPar.steps = xmlPar.stoutPar.steps;
    stPar.orthogDim = xmlPar.stoutPar.orthog;
    stPar.rho = xmlPar.stoutPar.rho;
    application.createModule<MGauge::StoutSmearing>("Vgauge_orig",stPar);

    // translat Vgauge field
    tPar.q = "Vgauge_orig";
    tPar.shift = xmlPar.configPar.source_loc;
    application.createModule<MUtilities::CShiftGaugeField>("Vgauge",tPar);

    // single prec cast
    if (xmlPar.mdwfPars.mixed) {
        MUtilities::GaugeSinglePrecisionCast::Par fPar;
        fPar.field = "Vgauge";
        application.createModule<MUtilities::GaugeSinglePrecisionCast>("VgaugeF",fPar);
    }

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

    std::vector<std::string> step0 = {};

    std::vector<std::string> qs = {};
    for (int k = 0; k < srcs.size(); k++) {
        int j = 0;
        if (k > 0) j = 1;

        for (qdx = 0; qdx < xmlPar.mdwfPars.masses.size(); qdx++) {
            std ::string gsrc, name_q1;

            if (xmlPar.smearPar.smear == 1) {
                gsrc = make_SmrSrc(application, xmlPar.smearPar, srcs[k], "trans_gauge", "_"+sources[k]);
                name_q1 = make_moebiusdwfpropagator(application, xmlPar.mdwfPars, qdx, gsrc, j, boundary, twist, "Vgauge", "VgaugeF", "_"+sources[k]);
            } else {
                name_q1 = make_moebiusdwfpropagator(application, xmlPar.mdwfPars, qdx, srcs[k], j, boundary, twist, "Vgauge", "VgaugeF", "_"+sources[k]);
            }
            qs.push_back(name_q1);

            MContraction::WardIdentity::Par WPar;
            WPar.prop = name_q1+"_5d";
            WPar.action = "MDWF_action_" + xmlPar.mdwfPar1.quark[qdx] + "_"+sources[0];
            WPar.mass = xmlPar.mdwfPar1.mass[qdx];
            if (xmlPar.smearPar.smear == 1) {
                WPar.source = gsrc;
            } else {
                WPar.source = srcs[k];
            }
            application.createModule<MContraction::WardIdentity>("WardId_"+xmlPar.mdwfPar1.quark[qdx]+sources[k],WPar);
            step0.push_back("WardId_"+xmlPar.mdwfPar1.quark[qdx]+sources[k]);
        }
    }

    // gradient flow
    MGradientFlow::WilsonFlow::Par gfPar;
    gfPar.gauge = "trans_gauge";
    gfPar.steps = 0; 
    gfPar.step_size = xmlPar.flowPar.step_size;
    gfPar.meas_interval = 1; 

    MGradientFlow::WilsonFermionFlow::Par wfPar;
        wfPar.gauge = "trans_gauge";
        wfPar.steps = 1; 
        wfPar.step_size = xmlPar.flowPar.step_size;
        wfPar.meas_interval = xmlPar.flowPar.meas_interval;
        wfPar.props = qs;
        wfPar.bc = xmlPar.flowPar.bc;
    if (xmlPar.flowPar.stout == 1) 
    {
        wfPar2.stoutSteps = std::to_string(xmlPar.stoutPar.steps);
        wfPar2.stoutRho = std::to_string(xmlPar.stoutPar.rho);
    }

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
    for (int k = 0; k < srcs.size(); k++) {
        for (qdx = 0; qdx < xmlPar.mdwfPars.masses.size(); qdx++) {
            for (rdx = 0; rdx < xmlPar.mdwfPars.masses.size(); rdx++) {
                std::string q1 = xmlPar.mdwfPars.quark[qdx], q2 = xmlPar.mdwfPars.quark[rdx];
                std::string name_q1 = q1+"_"+sources[k], name_q2 = q2+"_"+sources[k];
                step0.push_back(make_contraction(application, q1+q2+"_t"+tm+"_0.00", name_q1, name_q2, gamma_pairs, snk_str, "", true, false));
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
    const int base_meas_interval = (xmlPar.flowPar.meas_interval > 0) ? xmlPar.flowPar.meas_interval : 1;
    const int coarse_meas_interval = 4 * base_meas_interval;
    int coarsen_step = 1;
    if (xmlPar.flowPar.step_size > 0.0) {
        coarsen_step = static_cast<int>(std::ceil(xmlPar.flowPar.Tcoarsen / xmlPar.flowPar.step_size - 1e-12));
        if (coarsen_step < 1) coarsen_step = 1;
    }
    int last_contraction_step = 0;
    std::vector<std::string> outNames = {};

    for (int t = 1; t <= xmlPar.flowPar.steps; t++) {
        double ti = xmlPar.flowPar.step_size * t;
        int meas_interval = (t >= coarsen_step) ? coarse_meas_interval : base_meas_interval;

        // pass evolved propagators to the next step in the gradient flow
        if (t != 1) {
            wfPar.gauge = WFF + "_U";
            wfPar.props = outNames;
        }

        std::stringstream st; st << std::fixed << std::setprecision(2) << ti;
        fM1 = "FlowTime_t" + st.str(); WF1 = "WilsonFlow_t" + st.str(); 
        std::vector<std::string> stepi = {};

        std::string time = "_t" + st.str();

        outNames.clear();
        for (std::string q : qs) {
            outNames.push_back(q+time);
        }
        wfPar.outProps = outNames;

        // only perform contractions for set intervals along the evolution
        // always perform contractions at final flow time
        if ((t - last_contraction_step >= meas_interval) || (t == xmlPar.flowPar.steps)) {
            last_contraction_step = t;
            // add gauge observables important along gradient flow
            application.createModule<MGradientFlow::WilsonFermionFlow>(WF1,wfPar);
            
            // 2pt contractions
            for (int k = 0; k < srcs.size(); k++) {
                for (qdx = 0; qdx < xmlPar.mdwfPars.masses.size(); qdx++) {
                    for (rdx = 0; rdx < xmlPar.mdwfPars.masses.size(); rdx++) {
                        std::string q1 = xmlPar.mdwfPars.quark[qdx], q2 = xmlPar.mdwfPars.quark[rdx];
                        std::string flow_q1 = q1+"_"+sources[k]+time, flow_q2 = q2+"_"+sources[k]+time;
                        stepi.push_back(make_contraction(application, q1+q2+"_t"+tm+"_"+st.str(), flow_q1, flow_q2, gamma_pairs, snk_str, "", true, false));
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
    application.saveParameterFile(xmlPar.mesonPar.file_name+".xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
