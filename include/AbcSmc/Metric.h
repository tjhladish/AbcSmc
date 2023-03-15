
#ifndef ABCSMC_METRIC_H
#define ABCSMC_METRIC_H

namespace ABC {

    struct Metric {
        virtual std::string get_name() const = 0;
        virtual std::string get_short_name() const = 0;
        virtual bool is_integral() const = 0;
        virtual double get_obs_val() const = 0;        
    };

    template <typename NT> requires (std::is_integral_v<NT> or std::is_floating_point_v<NT>)
    class TMetric : public Metric {
        public:
            TMetric() {};
            TMetric(std::string s, std::string ss, double val) : name(s), short_name(ss), obs_val(val) {};

            std::string get_name() const override { return name; }
            std::string get_short_name() const override { if (short_name == "") { return name; } else { return short_name; } }
            bool is_integral() const override { if constexpr (std::integral<NT>) { return true; } else { return false; } }
            double get_obs_val() const override { return obs_val; }

        private:
            std::string name;
            std::string short_name;
            double obs_val;
    };
}

#endif // ABCSMC_METRIC_H