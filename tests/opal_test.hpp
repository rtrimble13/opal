// Minimal single-header test framework.
#pragma once

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace opal_test {

struct Registry {
    struct Case {
        std::string name;
        std::function<void()> fn;
    };
    static Registry& instance() {
        static Registry r;
        return r;
    }
    std::vector<Case> cases;
    int failures = 0;
    int checks = 0;
    std::string current;
};

struct Registrar {
    Registrar(const std::string& name, std::function<void()> fn) {
        Registry::instance().cases.push_back({name, std::move(fn)});
    }
};

inline void check_close(double actual, double expected, double tol,
                        const char* expr, const char* file, int line) {
    auto& r = Registry::instance();
    ++r.checks;
    bool ok = std::isfinite(actual) && std::fabs(actual - expected) <= tol;
    if (!ok) {
        ++r.failures;
        std::printf("FAIL [%s] %s:%d\n  %s\n  actual=%.10g expected=%.10g tol=%.3g diff=%.3g\n",
                    r.current.c_str(), file, line, expr, actual, expected, tol,
                    std::fabs(actual - expected));
    }
}

inline void check_true(bool cond, const char* expr, const char* file, int line) {
    auto& r = Registry::instance();
    ++r.checks;
    if (!cond) {
        ++r.failures;
        std::printf("FAIL [%s] %s:%d\n  %s\n", r.current.c_str(), file, line, expr);
    }
}

inline int run_all() {
    auto& r = Registry::instance();
    for (auto& c : r.cases) {
        r.current = c.name;
        std::printf("RUN  %s\n", c.name.c_str());
        try {
            c.fn();
        } catch (const std::exception& e) {
            ++r.failures;
            std::printf("FAIL [%s] unexpected exception: %s\n", c.name.c_str(),
                        e.what());
        }
    }
    std::printf("\n%d checks, %d failures\n", r.checks, r.failures);
    return r.failures == 0 ? 0 : 1;
}

}  // namespace opal_test

#define TEST_CASE(name)                                                       \
    static void test_fn_##name();                                             \
    static opal_test::Registrar reg_##name(#name, test_fn_##name);            \
    static void test_fn_##name()

#define CHECK_CLOSE(actual, expected, tol)                                    \
    opal_test::check_close((actual), (expected), (tol), #actual, __FILE__, __LINE__)

#define CHECK_TRUE(cond) opal_test::check_true((cond), #cond, __FILE__, __LINE__)

#define CHECK_THROWS(stmt)                                                    \
    do {                                                                      \
        bool threw_ = false;                                                  \
        try {                                                                 \
            stmt;                                                             \
        } catch (...) {                                                       \
            threw_ = true;                                                    \
        }                                                                     \
        opal_test::check_true(threw_, "expected throw: " #stmt, __FILE__,     \
                              __LINE__);                                      \
    } while (0)
