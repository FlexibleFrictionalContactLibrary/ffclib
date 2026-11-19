// cfc_solve: solves cfc/fclib files (FORMAT.md) with Clarabel, QOCO, ConicXX, SCS and cuOpt and writes
// one JSON line per run (schema: README.md, "Result schema"). The check of an export's residuals is
// cfc_info --check.
//
//   cfc_solve --list
//   cfc_solve [--solver a,b] [--form V,F] [--tol 1e-6,1e-8] [--repeat K] [--out results.jsonl] FILE...

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

#include <Eigen/Core>

#include <cfc/cone_program.hpp>
#include <cfc/file.hpp>
#include <cfc/json.hpp>
#include <cfc/residuals.hpp>

#include "solvers.hpp"

namespace {

using Obj = cfc::json::Object;

// ---------- platform ----------

std::string firstLine(const std::string& cmd) {
    std::string line;
    if (FILE* p = popen(cmd.c_str(), "r")) {
        char buf[512];
        if (std::fgets(buf, sizeof buf, p)) line = buf;
        pclose(p);
    }
    while (!line.empty() && (line.back() == '\n' || line.back() == ' ')) line.pop_back();
    return line;
}

Obj platform() {
    char host[256] = {0};
    gethostname(host, sizeof host - 1);
    std::string cpu;
    std::ifstream info("/proc/cpuinfo");
    for (std::string line; std::getline(info, line);)
        if (line.rfind("model name", 0) == 0) {
            cpu = line.substr(line.find(':') + 2);
            break;
        }
    const char* omp = std::getenv("OMP_NUM_THREADS");
    Obj o;
    o.str("host", host).str("cpu", cpu).integer("threads", omp ? std::atoi(omp) : 0);
    o.str("gpu", firstLine("nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null"));
    return o;
}

// ---------- options ----------

std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream in(s);
    for (std::string item; std::getline(in, item, ',');)
        if (!item.empty()) out.push_back(item);
    return out;
}

struct Options {
    bool list{false}, verbose{false};
    std::vector<std::string> solvers;
    std::vector<cfc::Form> forms{cfc::Form::V};
    std::vector<double> tols{1e-8};
    int max_iter{500};
    double time_limit{600};
    int repeat{3};
    bool zero_shift{false};
    std::string out;
    std::vector<std::string> files;
};

[[noreturn]] void usage(int code) {
    std::cout << "usage: cfc_solve [--list] [--solver a,b] [--form V,F] [--tol 1e-6,1e-8] [--max-iter N] [--time-limit S]\n"
                 "                 [--repeat K] [--shift engine|zero] [--out results.jsonl] [--verbose] FILE...\n"
                 "  --repeat K   K timed solves after one untimed warm-up (K = 1: no warm-up); medians are reported\n";
    std::exit(code);
}

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) usage(2);
            return argv[++i];
        };
        if (a == "--list")
            o.list = true;
        else if (a == "--verbose")
            o.verbose = true;
        else if (a == "--solver")
            o.solvers = split(value());
        else if (a == "--form") {
            o.forms.clear();
            for (const auto& f : split(value())) o.forms.push_back(f == "F" || f == "f" ? cfc::Form::F : cfc::Form::V);
        } else if (a == "--tol") {
            o.tols.clear();
            for (const auto& t : split(value())) o.tols.push_back(std::stod(t));
        } else if (a == "--max-iter")
            o.max_iter = std::stoi(value());
        else if (a == "--time-limit")
            o.time_limit = std::stod(value());
        else if (a == "--repeat")
            o.repeat = std::max(1, std::stoi(value()));
        else if (a == "--shift")
            o.zero_shift = value() == "zero";
        else if (a == "--out")
            o.out = value();
        else if (a == "-h" || a == "--help")
            usage(0);
        else if (!a.empty() && a[0] == '-')
            usage(2);
        else
            o.files.push_back(a);
    }
    if (o.solvers.empty()) o.solvers = cfcbench::available();
    return o;
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return n == 0 ? std::nan("") : (n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]));
}

std::string family(const cfc::Record& rec, const std::string& path) {
    if (rec.meta.present && !rec.meta.scene.empty()) return rec.meta.scene;
    return std::filesystem::path(path).parent_path().filename().string();
}

std::string problemName(const std::string& path) {
    std::string stem = std::filesystem::path(path).filename().string();
    for (const char* ext : {".cfc.h5", ".hdf5", ".h5"})
        if (stem.size() > std::strlen(ext) && stem.compare(stem.size() - std::strlen(ext), std::string::npos, ext) == 0) return stem.substr(0, stem.size() - std::strlen(ext));
    return stem;
}

Obj problemJson(Obj o, const cfc::Record& rec, const std::string& path) {
    const cfc::Problem& pb = rec.problem;
    o.str("problem", problemName(path)).str("file", path).str("family", family(rec, path)).integer("step", rec.meta.step);
    o.integer("nc", pb.nc()).integer("n", pb.n()).integer("m", 3 * pb.nc()).integer("p", pb.p());
    long long compliant = 0;
    for (Eigen::Index k = 0; k < pb.C.outerSize(); ++k)
        for (cfc::SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
            if (it.value() != 0) {
                ++compliant;
                break;
            }
    o.integer("n_compliant", compliant).integer("n_frictionless", (pb.mu.array() == 0).count());
    o.integer("nnz_H", pb.H.nonZeros()).integer("nnz_G", pb.G.nonZeros()).integer("nnz_C", pb.C.nonZeros());
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    Eigen::setNbThreads(1);
    const Options opt = parse(argc, argv);
    if (opt.list) {
        for (const auto& b : cfcbench::available()) std::cout << b << ' ' << cfcbench::version(b) << '\n';
        return 0;
    }
    if (opt.files.empty()) usage(2);

    std::ofstream file;
    if (!opt.out.empty()) file.open(opt.out, std::ios::app);
    std::ostream& out = opt.out.empty() ? std::cout : file;
    const Obj plat = platform();

    // cuOpt: CUDA context and cuDSS start-up once, outside all timings.
    double cudaInit = -1;
    std::string cudaReason;
    if (std::find(opt.solvers.begin(), opt.solvers.end(), "cuopt") != opt.solvers.end()) {
        cudaInit = cfcbench::cudaWarmup(&cudaReason);
        if (cudaInit < 0)
            std::cerr << "cuopt: no usable GPU (" << cudaReason << ")\n";
        else
            std::cerr << "cuopt: CUDA/cuDSS start-up " << cudaInit << " s (excluded from timings)\n";
    }

    for (const auto& path : opt.files) {
        const auto r0 = std::chrono::steady_clock::now();
        const cfc::Record rec = cfc::read(path);
        const double tRead = std::chrono::duration<double>(std::chrono::steady_clock::now() - r0).count();
        const cfc::Problem& pb = rec.problem;
        const cfc::Vector shift = opt.zero_shift ? cfc::Vector::Zero(pb.nc()) : pb.shift;
        cfc::json::Object git;
        git.str("bench", CFCBENCH_GIT).str("producer", rec.meta.producer_version);

        for (const cfc::Form form : opt.forms) {
            const auto b0 = std::chrono::steady_clock::now();
            const cfc::ConeProgram cp = cfc::buildConeProgram(pb, shift, form);
            const double tBuild = std::chrono::duration<double>(std::chrono::steady_clock::now() - b0).count();
            for (const auto& solver : opt.solvers) {
                if (!cfcbench::acceptsForm(solver, form)) continue;
                for (const double tol : opt.tols) {
                    cfcbench::Settings st;
                    st.tol = tol;
                    st.max_iter = opt.max_iter;
                    st.time_limit = opt.time_limit;
                    st.verbose = opt.verbose;
                    std::vector<double> setup, solve, reported;
                    cfcbench::Result res;
                    const int runs = opt.repeat > 1 ? opt.repeat + 1 : 1;
                    for (int k = 0; k < runs; ++k) {
                        res = cfcbench::solve(solver, cp, st);
                        if (runs > 1 && k == 0) continue;  // warm-up
                        setup.push_back(res.t_setup);
                        solve.push_back(res.t_solve);
                        reported.push_back(res.t_solver);
                        if (!res.ok) break;  // do not repeat failures
                    }
                    Obj o = problemJson(Obj(), rec, path);
                    o.str("form", cfc::toString(form)).str("shift", opt.zero_shift ? "zero" : "engine").str("solver", solver).str("solver_version", cfcbench::version(solver));
                    o.obj("platform", plat);
                    Obj settings;
                    settings.num("tol_target", tol).integer("max_iter", opt.max_iter).num("time_limit", opt.time_limit).integer("repeat", static_cast<long long>(solve.size()));
                    o.obj("settings", settings);
                    o.str("status", res.status).boolean("ok", res.ok).boolean("solved", res.solved).integer("iterations", res.iterations);
                    o.num("t_read", tRead)
                        .num("t_build", tBuild)
                        .num("t_setup", median(setup))
                        .num("t_solve", median(solve))
                        .num("t_total", median(setup) + median(solve))
                        .num("t_solver", median(reported));
                    o.nums("t_solve_samples", solve).num("t_cuda_init", solver == "cuopt" ? cudaInit : std::nan(""));
                    if (res.ok && res.x.size() == cp.P.cols() && res.z.size() == cp.A.rows()) {
                        const cfc::Solution sol = cfc::recoverSolution(pb, form, res.x, res.z);
                        o.obj("residuals", cfc::json::residuals(cfc::residuals(pb, sol, shift)));
                    } else {
                        o.raw("residuals", "null");
                    }
                    o.num("objective", res.objective).str("message", res.message).obj("git", git);
                    out << o.dump() << std::endl;
                    if (opt.verbose || opt.out.size())
                        std::cerr << problemName(path) << " " << cfc::toString(form) << " " << solver << " tol " << tol << ": " << res.status << ", " << res.iterations << " it, "
                                  << median(setup) + median(solve) << " s\n";
                }
            }
        }
    }
    return 0;
}
