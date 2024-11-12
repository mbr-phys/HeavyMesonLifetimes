#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

#include "Inputs.hpp"
#include "Functions.hpp"
//#include "RHQ.hpp"

//#include <sstream>

using namespace Grid;
using namespace Hadrons;

// set fermion boundary conditions to be periodic space, antiperiodic time.
// these will likely need to be generalised, e.g. in xml, at some point
std::string boundary = "1 1 1 -1";
std::string twist = "0. 0. 0. 0.";

struct TestPar
{
    TestInputs::RunPar      runPar;
    TestInputs::ConfigPar   configPar;
    TestInputs::MdwfPar     mdwfParh;
    TestInputs::DwfPar      dwfPars;
    TestInputs::MesonPar    mesonPar;
    TestInputs::FlowPar     flowPar;
    TestInputs::StoutPar    stoutPar;
    TestInputs::SmearPar    smearPar;
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

    TestPar testPar;
    XmlReader reader(parFilename);

    read(reader, "runPar",      testPar.runPar);
    read(reader, "configPar",   testPar.configPar);
    read(reader, "hdwfPar",     testPar.mdwfParh);
    read(reader, "sdwfPar",     testPar.dwfPars);
    read(reader, "mesonPar",    testPar.mesonPar);
    read(reader, "flowPar",     testPar.flowPar);
    read(reader, "stoutPar",    testPar.stoutPar);
    read(reader, "smearPar",    testPar.smearPar);

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
    globalPar.trajCounter.start        = testPar.configPar.begin;
    globalPar.trajCounter.end          = testPar.configPar.end;
    globalPar.trajCounter.step         = testPar.configPar.step;
    globalPar.runId                    = testPar.runPar.runId;
    globalPar.database.applicationDb   = testPar.runPar.appDb;
    globalPar.database.restoreSchedule = testPar.runPar.restore;
    if (testPar.runPar.scheduler == "naive") {
        globalPar.scheduler.schedulerType = "naive";}
    application.setPar(globalPar);
    
    // gauge field
    MIO::LoadNersc::Par gauge;
    gauge.file = testPar.configPar.fileStem;
    application.createModule<MIO::LoadNersc>("gauge",gauge);

    // translat gauge field
    MGauge::Translat::Par tPar;
    tPar.gauge = "gauge";
    tPar.xvec = coordVec(testPar.configPar.source_loc);
    application.createModule<MGauge::Translat>("trans_gauge",tPar);

    // stout smear
    MGauge::StoutSmearing::Par stPar;
    stPar.gauge = "trans_gauge";
    stPar.steps = testPar.stoutPar.steps;
    stPar.orthogDim = testPar.stoutPar.orthog;
    stPar.rho = testPar.stoutPar.rho;
    application.createModule<MGauge::StoutSmearing>("Vgauge",stPar);

    // single prec cast
    MUtilities::GaugeSinglePrecisionCast::Par fPar;
    fPar.field = "trans_gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("trans_gaugeF",fPar);

    std::string gamma_pairs = "";
    if (testPar.mesonPar.gamma_snk == "") {
        gamma_pairs = "all";}
    else {
        std::vector<std::string> gamma_src_vec = split_gammas(testPar.mesonPar.gamma_src,' ');
        std::vector<std::string> gamma_snk_vec = split_gammas(testPar.mesonPar.gamma_snk,' ');
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
    for (int time : testPar.configPar.sources) {
        std::stringstream ss; ss << time;
        testPar.configPar.source_loc = "0 0 0 " + ss.str();
        sources.push_back(ss.str());
        srcs.push_back(make_Z2src(application, time, ""));
    }

    std::vector<std::string> hqs = {};
    std::vector<std::string> sqs = {};
    std::vector<std::string> step0 = {};
    for (int k = 0; k < srcs.size(); k++) {
        int j = 0;
        if (k > 0) j = 1;
        // heavy quark
        std::string name_hq1 = make_moebiusdwfpropagator(application, testPar.mdwfParh, srcs[k], j, boundary, twist, "Vgauge", "_"+sources[k]);
        hqs.push_back(name_hq1);

        MContraction::WardIdentity::Par WPar;
        WPar.prop = name_hq1+"_5d";
        WPar.action = "MDWF_action_" + testPar.mdwfParh.quark + "_"+sources[0];
        WPar.mass = testPar.mdwfParh.mass;
        WPar.source = srcs[k];
        application.createModule<MContraction::WardIdentity>("WardId_"+testPar.mdwfParh.quark+sources[k],WPar);
        step0.push_back("WardId_"+testPar.mdwfParh.quark+sources[k]);

        std ::string gsrc = make_SmrSrc(application, testPar.smearPar, srcs[k], "trans_gauge", "_"+sources[k]);
        // strange quark
        //std::string name_sq1 = make_mixedmoebiusdwfpropagator(application, testPar.mdwfPars, gsrc, j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+sources[k]);
        std::string name_sq1 = make_mixeddwfpropagator(application, testPar.dwfPars, gsrc, j, boundary, twist, "trans_gauge", "trans_gaugeF", "_"+sources[k]);
        sqs.push_back(name_sq1);

        MContraction::WardIdentity::Par SPar;
        SPar.prop = name_sq1+"_5d";
        SPar.action = "DWF_action_" + testPar.dwfPars.quark + "_" + sources[0];
        SPar.mass = testPar.dwfPars.mass;
        SPar.source = srcs[k];
        application.createModule<MContraction::WardIdentity>("WardId_"+testPar.dwfPars.quark+sources[k],SPar);
        step0.push_back("WardId_"+testPar.dwfPars.quark+sources[k]);
    }
    std::vector<std::string> qs = hqs; 
    qs.insert(qs.end(), sqs.begin(), sqs.end());

    std::string fileSt = fileStem(testPar.mesonPar.cont_name);

    std::string snk_mom = "0 0 0"; 
    std::string snk_str = make_PtSnk(application, snk_mom, "");

    int Tmax = testPar.configPar.Tmax;

    std::string flowMod = "FlowTime_t";
    std::string WFF = "WilsonFlow_t";
    std::string fM1(flowMod), WF1(WFF);
    for (int q = 0; q < hqs.size(); q++) {
        std::string name_sq1 = sqs[q], name_hq1 = hqs[q], tm = sources[q];
        step0.push_back(make_contraction(application, "sh_SP_t"+tm+"_0.00", name_sq1, name_hq1, gamma_pairs, snk_str, "", true, false));
        step0.push_back(make_contraction(application, "ss_SP_t"+tm+"_0.00", name_sq1, name_sq1, gamma_pairs, snk_str, "", true, false));
        step0.push_back(make_contraction(application, "hh_PP_t"+tm+"_0.00", name_hq1, name_hq1, gamma_pairs, snk_str, "", true, false));
        if (testPar.mesonPar.ss == 1) {
            std::string name_sqsm = make_SmrSrc(application, testPar.smearPar, name_sq1, "trans_gauge", "_smr");
            step0.push_back(make_contraction(application, "sh_SS_t"+tm+"_0.00", name_sqsm, name_hq1, gamma_pairs, snk_str, "", true, false));
            step0.push_back(make_contraction(application, "ss_SS_t"+tm+"_0.00", name_sqsm, name_sqsm, gamma_pairs, snk_str, "", true, false));
        }
    }
    std::vector<std::string> corrgroups = {make_corrgroup(application, flowMod+"0.00",step0)};

    MIO::WriteCorrelatorGroup::Par WCPar;
    WCPar.contractions = corrgroups;
    WCPar.output = fileSt + "LumiCheck";
    application.createModule<MIO::WriteCorrelatorGroup>("WriteToFile",WCPar);

    // execution
    application.saveParameterFile("test.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
