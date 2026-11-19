#pragma once

#include <string>

#include "problem.hpp"

// cfc files (FORMAT.md): HDF5 files in the layout of fclib's global problem (/fclib_global, /solution,
// /guesses) with two optional extensions that fclib's readers ignore:
//
//   /fclib_global/C                       compliance C (p x p, fclib's matrix layout)
//   /fclib_global/vectors/desaxce_shift   frozen shift of the convexified law (nc)
//   /cfc                                  metadata (Metadata)
//
// Every plain fclib file is a cfc file (C = 0, shift = 0, no metadata).
namespace cfc {

struct WriteOptions {
    int compression{4};  // gzip level of datasets with at least 1024 values (with shuffle); 0: none
};

// Writes `rec` to `path`, replacing an existing file. fclib needs at least one contact. Throws
// std::runtime_error on failure, std::invalid_argument for an inconsistent record.
void write(const std::string& path, const Record& rec, const WriteOptions& options = {});

// Reads a cfc file or a plain fclib global problem (3D). Throws std::runtime_error if the file holds no
// global problem or has an unsupported /cfc version.
Record read(const std::string& path);

}  // namespace cfc
