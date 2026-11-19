#pragma once

#include <optional>
#include <string>
#include <vector>

#include "problem.hpp"

// Files in fclib's HDF5 layout of the global problem (fclib.h of FrictionalContactLibrary/fclib) with one
// extension, the compliance /fclib_global/C (FORMAT.md). Every 3D global problem written by fclib is read,
// with C = 0; other groups in a file are ignored.
namespace ffclib {

// fclib's free-text description of a problem (/fclib_global/info).
struct Info {
    std::string title;
    std::string description;
    std::string math_info;  // known mathematical properties; producers state the compliance here
};

// Everything a file holds.
struct Record {
    Problem problem;
    Info info;
    std::optional<Solution> solution;  // /solution, e.g. the producer's
    std::vector<Solution> guesses;     // /guesses, e.g. the producer's warm start
};

// Throws std::runtime_error if the file cannot be read or holds no valid 3D global problem (check()).
Record read(const std::string& path);

// Writes `record` to `path`, replacing an existing file. Throws std::invalid_argument for an invalid record
// (fclib needs at least one contact), std::runtime_error if writing fails.
void write(const std::string& path, const Record& record);

}  // namespace ffclib
