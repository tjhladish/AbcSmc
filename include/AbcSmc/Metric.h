
#ifndef ABCSMC_METRIC_H
#define ABCSMC_METRIC_H

#include <type_traits>
#include <memory>

namespace ABC {

    struct Metric {
        Metric(std::string s, std::string ss) : name(s), short_name(ss) {}
        std::string get_name() const { return name; }
        std::string get_short_name() const { if (short_name == "") { return name; } else { return short_name; } }
        virtual bool is_integral() const = 0;
        virtual double get_obs_val() const = 0;
        
        private:
            std::string name;
            std::string short_name;
    };

    typedef std::shared_ptr<const Metric> MetricPtr;
    typedef std::vector<MetricPtr> MetricVec;

    // Type'd metric (as in, integer or float typed)
    template <typename NT>
    class TMetric : public Metric {
        public:
            TMetric(std::string s, std::string ss, double val) : Metric(s, ss), obs_val(val) {};

            
            bool is_integral() const override { if constexpr (std::is_integral_v<NT>) { return true; } else { return false; } }
            double get_obs_val() const override { return obs_val; }

        private:
            std::string name;
            std::string short_name;
            double obs_val;
    };
}

#endif // ABCSMC_METRIC_H