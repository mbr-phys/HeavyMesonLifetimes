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
                                        bool, restore);
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
                                        int, Tcoarsen,
                                        std::string, maxTau,
                                        int, bc); 
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

