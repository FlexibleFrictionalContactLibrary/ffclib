#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "residuals.hpp"

// Minimal JSON output (one object per line) for tools that report problems and residuals.
namespace cfc::json {

inline std::string quote(const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\')
            out += '\\', out += c;
        else if (c == '\n')
            out += "\\n";
        else if (static_cast<unsigned char>(c) < 0x20)
            out += ' ';
        else
            out += c;
    }
    return out + "\"";
}

// Shortest round-trip representation; NaN and infinities as null.
inline std::string number(double x) {
    if (!std::isfinite(x)) return "null";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.17g", x);
    return buf;
}

// Object with the keys in insertion order.
class Object {
   public:
    Object& raw(const std::string& k, const std::string& v) {
        m_items.emplace_back(k, v);
        return *this;
    }
    Object& str(const std::string& k, const std::string& v) { return raw(k, quote(v)); }
    Object& num(const std::string& k, double v) { return raw(k, number(v)); }
    Object& integer(const std::string& k, long long v) { return raw(k, std::to_string(v)); }
    Object& boolean(const std::string& k, bool v) { return raw(k, v ? "true" : "false"); }
    Object& obj(const std::string& k, const Object& o) { return raw(k, o.dump()); }
    Object& nums(const std::string& k, const std::vector<double>& v) {
        std::string s = "[";
        for (std::size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + number(v[i]);
        return raw(k, s + "]");
    }
    std::string dump() const {
        std::string s = "{";
        for (std::size_t i = 0; i < m_items.size(); ++i) s += (i ? ", " : "") + quote(m_items[i].first) + ": " + m_items[i].second;
        return s + "}";
    }

   private:
    std::vector<std::pair<std::string, std::string>> m_items;
};

// The residuals with rel_max = convexifiedMax() / velocity_scale.
inline Object residuals(const Residuals& r) {
    Object o;
    o.num("balance", r.balance).num("bilateral", r.bilateral).num("natural_map", r.natural_map).num("natural_map_coulomb", r.natural_map_coulomb);
    o.num("cone_primal", r.cone_primal).num("cone_dual", r.cone_dual).num("compl", r.complementarity).num("velocity_scale", r.velocity_scale);
    o.num("rel_max", r.convexifiedMax() / r.velocity_scale);
    return o;
}

}  // namespace cfc::json
