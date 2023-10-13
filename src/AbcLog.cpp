#include <AbcSmc/AbcLog.h>

// TODO: longer term goal is to have these *not* have an AbcSmc argument
// and instead have the specific elements they need passed in as arguments

void AbcLog::print_stats(
    const std::string &str1, const std::string &str2,
    const double val1, const double val2, const double delta, const double pct_chg,
    const std::string &tail, ostream &os
) {
    os << "    " + str1 + ", " + str2 + "  ( delta, % ): "  << setw(WIDTH) << val1 << ", " << setw(WIDTH) << val2
                                                      << " ( " << setw(WIDTH) << delta << ", " << setw(WIDTH) << pct_chg  << "% )\n" + tail;
}


void AbcLog::_print_particle_table_header(
    AbcSmc * abc,
    std::ostream &os
) {
    for (size_t i = 0; i < abc->npar(); i++) { os << setw(WIDTH) << abc->_model_pars[i]->get_short_name(); } os << " | ";
    for (size_t i = 0; i < abc->nmet(); i++) { os << setw(WIDTH) << abc->_model_mets[i]->get_short_name(); } os << std::endl;
}

void AbcLog::report_convergence_data(
    AbcSmc *abc,
    const size_t set_t,
    std::ostream &os
) {
    if( abc->_predictive_prior.size() <= set_t ) {
        os << "ERROR: attempting to report stats for set " << set_t << ", but data aren't available. " << std::endl
             << "       This can happen if --process is called on a database that is not ready to be processed." << std::endl;
        exit(-214);
    }

    Mat2D par_values = abc->_particle_parameters[set_t](abc->_predictive_prior[set_t], Eigen::placeholders::all);
    Row current_means = par_values.colwise().mean();
    Row last_means;
    if (set_t > 0) { // if there was a previous set, also extract last_means
        Mat2D last_par_values = abc->_particle_parameters[set_t-1](abc->_predictive_prior[set_t-1], Eigen::placeholders::all);
        last_means = last_par_values.colwise().mean();
    }

    os << double_bar << std::endl;
    if (set_t == 0) {
        os << "Predictive prior summary statistics:\n";
    } else {
        os << "Convergence data for predictive priors:\n";
    }
    for (size_t parIdx = 0; parIdx < abc->_model_pars.size(); parIdx++) {
        const ABC::Parameter* par = abc->_model_pars[parIdx];
        const double current_stdev = sqrt(abc->_doubled_variance[set_t][parIdx]/2.0);
        const double prior_mean = par->get_mean();
        const double prior_mean_delta = current_means[parIdx] - prior_mean;
        const double prior_mean_pct_chg = prior_mean != 0 ? 100 * prior_mean_delta / prior_mean : INFINITY;

        const double prior_stdev = par->get_sd();
        const double prior_stdev_delta = current_stdev - prior_stdev;
        const double prior_stdev_pct_chg = prior_stdev != 0 ? 100 * prior_stdev_delta / prior_stdev : INFINITY;
        os << "  Par " << parIdx << ": \"" << par->get_name() << "\"\n";
        os << "  Means:\n";
        print_stats("Prior", "current", prior_mean, current_means[parIdx], prior_mean_delta, prior_mean_pct_chg, "", os);

        if (set_t != 0) {
            double delta = current_means[parIdx] - last_means[parIdx];
            double pct_chg = last_means[parIdx] != 0 ? 100 * delta / last_means[parIdx] : INFINITY;
            print_stats("Last", " current", last_means[parIdx], current_means[parIdx], delta, pct_chg, "\n", os);
        }

        os << "  Standard deviations:\n";
        print_stats("Prior", "current", prior_stdev, current_stdev, prior_stdev_delta, prior_stdev_pct_chg, "\n", os);

        if (set_t != 0) {
            double last_stdev = sqrt(abc->_doubled_variance[set_t-1][parIdx]/2.0);
            double delta = current_stdev - last_stdev;
            double pct_chg = last_stdev != 0 ? 100 * delta / last_stdev : INFINITY;
            print_stats("Last", " current", last_stdev, current_stdev, delta, pct_chg, "\n", os);
        }
    }
}

void AbcLog::filtering_report(
    AbcSmc * abc,
    const size_t t,
    const Mat2D &posterior_pars, // rows = samples, by rank; cols = parameters, by order
    const Mat2D &posterior_mets, // rows = samples, by rank; cols = metrics, by order
    std::ostream &os
) {
    os << double_bar << std::endl << "Set " << t << std::endl << double_bar << std::endl;

    os << "Observed:" << std::endl;
    AbcLog::_print_particle_table_header(abc, os);
    for (size_t i = 0; i < static_cast<size_t>(posterior_pars.cols()); i++) { os << setw(WIDTH) << "---"; } os << " | ";
    for (auto modmet : abc->_model_mets) { os << setw(WIDTH) << modmet->get_obs_val(); } os << std::endl;

    os << "Normalized RMSE for metric means (lower is better):  " << ABC::calculate_nrmse(posterior_mets, abc->_met_vals) << std::endl;
    os << "Posterior means:" << std::endl;
    AbcLog::_print_particle_table_header(abc, os);
    for (auto meanpar : posterior_pars.colwise().mean()) { os << setw(WIDTH) << meanpar; } os << " | ";
    for (auto meanmet : posterior_mets.colwise().mean()) { os << setw(WIDTH) << meanmet; } os << std::endl;

    os << "Posterior medians:" << std::endl;
    AbcLog::_print_particle_table_header(abc, os);
    for (auto parcol : posterior_pars.colwise()) { os << setw(WIDTH) << ABC::median(parcol); } os << " | ";
    for (auto metcol : posterior_mets.colwise()) { os << setw(WIDTH) << ABC::median(metcol); } os << std::endl;

    os << "Best five:" << std::endl;
    AbcLog::_print_particle_table_header(abc, os);
    auto partop = posterior_pars.topRows(5);
    auto mettop = posterior_mets.topRows(5);
    for (size_t q = 0; q < 5; q++) {
        for (auto par : partop.row(q)) { os << setw(WIDTH) << par; } os << " | ";
        for (auto met : mettop.row(q)) { os << setw(WIDTH) << met; } os << std::endl;
    }

    os << "Worst five:" << std::endl;
    AbcLog::_print_particle_table_header(abc, os);
    auto parbot = posterior_pars.bottomRows(5);
    auto metbot = posterior_mets.bottomRows(5);
    for (size_t q = 0; q < 5; q++) {
        for (auto par : parbot.row(q)) { os << setw(WIDTH) << par; } os << " | ";
        for (auto met : metbot.row(q)) { os << setw(WIDTH) << met; } os << std::endl;
    }

}