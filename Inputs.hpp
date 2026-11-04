#pragma once

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

namespace XmlInputs
{
    class RunPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(RunPar, 
                                        std::string, runId,
                                        std::string, appDb,
                                        std::string, scheduler,
                                        bool, restore,
                                        std::string, flowChaining,  // "step" or "interval"
                                        unsigned int, eyes,         // 1: dim-6 DeltaB=0 eyes and condensates (stochasticPar) in this run
                                        unsigned int, saveWalls,    // 1: save the unflowed wall propagators for FermionFlow_Eyes
                                        std::string, wallStem,      // directory of the saved wall propagators (<wallStem>/<prop>.<traj>.bin)
                                        unsigned int, wallIo32);    // 1: save/load the wall propagators in single precision
    };

    class ConfigPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ConfigPar,
                                        std::string,  fileStem,
                                        unsigned int, begin,
                                        unsigned int, end,
                                        unsigned int, step,
                                        std::string, source_loc,
                                        std::vector<int>, sources,
                                        std::vector<int>, deltaTs,
                                        unsigned int, Tmax);
    };
    
    class DwfPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(DwfPar,
                                        std::string, quark,
                                        double, mass,
                                        double, M5,
                                        unsigned int, Ls,
                                        double, residual,
                                        unsigned int, maxiter,
                                        unsigned int, maxiterOut,
                                        unsigned int, maxiterIn,
                                        std::string, prop);
    };

    class MdwfPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MdwfPar,
                                        std::string, quark,
                                        double, mass,
                                        double, M5,
                                        double, b,
                                        double, c,
                                        unsigned int, Ls,
                                        double, residual,
                                        unsigned int, maxiter,
                                        unsigned int, maxiterOut,
                                        unsigned int, maxiterIn,
                                        std::string, prop);
    };

    class MesonPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                        std::string, file_name,
                                        std::string, gamma_src,
                                        std::string, gamma_snk,
                                        std::string, rhq_impr,
                                        unsigned int, light,
                                        unsigned int, n2_min,
                                        unsigned int, n2_max,
                                        unsigned int, ss,
                                        std::string, moms);
    };

    class SmearPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(SmearPar,
                                        unsigned int, smear,
                                        double, width,
                                        int, iterations,
                                        int, orthog);
    };

    class FlowPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FlowPar,
                                        std::string, action,
                                        std::string, output,
                                        int, steps,
                                        double, step_size,
                                        int, meas_interval,
                                        double, Tcoarsen,
                                        std::string, maxTau,
                                        int, bc); 
    };

    class HQEPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(HQEPar,
                                        std::string, probeGammas,   // source/sink of all 3pts: Gamma5 and/or GammaTGamma5
                                        std::string, vectorGammas,  // B* 3pts (sink = source): "none" or a subset of GammaX GammaY GammaZ
                                        bool, twoQuark,             // HeavyMeson3pt two-quark (HQE) instances
                                        std::string, vertices,      // two-quark vertex sets, empty = "all ops"
                                        bool, dim7,                 // WeakNonEye3ptDeriv dim-7 DeltaB=2 instances (Lenz-Nierste)
                                        std::string, gammaPairs,    // dim-7 four-quark gamma pairs, empty = "diag"
                                        bool, atTau0);              // also run the HeavyMeson3pt/WeakNonEye3ptDeriv instances at zero flow time
    };

    class StochasticPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticPar,
                                        unsigned int, nSrcs,        // noise hits
                                        unsigned int, nsparse,      // sparse dilution per direction (nSrcs*nsparse^4 components)
                                        double, c_fl,               // flow-time improvement coefficient of the scalar condensate (as in SParameterGF)
                                        std::string, eyeGammaPairs); // WeakEye3pt gamma pairs: diag, all or "(G1 G2)..."; all for the penguin-type operators
    };

    class StoutPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(StoutPar,
                                        unsigned int, steps, 
                                        std::string, orthog,
                                        double, rho);
    };

}

