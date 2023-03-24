
#ifndef ABC_LOG_H
#define ABC_LOG_H

#include <iostream>

#include <AbcSmc/AbcSmc.h>

struct AbcLog {

    static void report_convergence_data(
        AbcSmc * abc, const size_t t,
        std::ostream & os = std::cerr
    );
    static void _print_particle_table_header(
        AbcSmc * abc,
        std::ostream & os = std::cerr
    );

    static void print_stats(
        const std::string & str1, const std::string & str2,
        const double val1, const double val2,
        const double delta, const double pct_chg,
        const string & tail,
        ostream& os = std::cerr
    );

    static void filtering_report(
        AbcSmc * abc,
        const size_t t,
        const Mat2D & posterior_pars, // rows = samples, by rank; cols = parameters, by order
        const Mat2D & posterior_mets, // rows = samples, by rank; cols = metrics, by order
        std::ostream & os = std::cerr
    );


    inline static const int WIDTH = 12;
    inline static const string double_bar = "=========================================================================================";

    private:
        AbcLog() {};

};

#endif