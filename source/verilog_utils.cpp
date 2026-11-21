#include <algorithm>
#include <multiplierless/verilog_utils.hpp>
#include <regex>
#include <set>
#include <stdexcept>
#include <utility>

namespace {

    auto split_lines(const std::string& text) -> std::vector<std::string> {
        std::vector<std::string> lines;
        size_t pos = 0;
        while (pos < text.size()) {
            auto const eol = text.find('\n', pos);
            lines.push_back(text.substr(pos, eol - pos));
            pos = (eol == std::string::npos) ? text.size() : eol + 1;
        }
        if (!text.empty() && text.back() == '\n') {
            lines.emplace_back();
        }
        return lines;
    }

    auto trimmed(const std::string& line) -> std::string {
        auto const start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos) return {};
        auto const end = line.find_last_not_of(" \t\r");
        return line.substr(start, end - start + 1);
    }

    auto code_part(const std::string& line) -> std::string {
        auto const comment = line.find("//");
        auto code = line.substr(0, comment);
        while (!code.empty()
               && (code.back() == ' ' || code.back() == '\t' || code.back() == '\r')) {
            code.pop_back();
        }
        return code;
    }

}  // namespace

auto csd_coeffs_to_fixed(const std::vector<std::string>& csd_list) -> FixedCoefficients {
    if (csd_list.empty()) {
        throw std::invalid_argument("at least one coefficient is required");
    }
    std::vector<std::pair<std::string, std::string>> split;
    size_t int_width = 0;
    size_t frac_width = 0;
    for (auto const& csd : csd_list) {
        auto const dot = csd.find('.');
        auto int_part = (dot == std::string::npos) ? csd : csd.substr(0, dot);
        auto frac_part = (dot == std::string::npos) ? std::string{} : csd.substr(dot + 1);
        if (int_part.empty()) {
            throw std::invalid_argument("CSD string must have an integer digit");
        }
        int_width = std::max(int_width, int_part.size());
        frac_width = std::max(frac_width, frac_part.size());
        split.emplace_back(std::move(int_part), std::move(frac_part));
    }
    std::vector<std::string> digits;
    digits.reserve(split.size());
    for (auto const& [int_part, frac_part] : split) {
        std::string raw(int_width - int_part.size(), '0');
        raw += int_part;
        raw += frac_part;
        raw.append(frac_width - frac_part.size(), '0');
        digits.push_back(std::move(raw));
    }
    return {std::move(digits), static_cast<int>(int_width + frac_width) - 1,
            static_cast<int>(frac_width)};
}

auto fix_direct_verilog(std::string verilog, const std::vector<std::string>& port_names)
    -> std::string {
    auto lines = split_lines(verilog);
    std::vector<size_t> port_lines;
    bool in_ports = false;
    for (size_t i = 0; i < lines.size(); ++i) {
        auto const line = trimmed(lines[i]);
        if (!in_ports && line.starts_with("module ") && lines[i].find('(') != std::string::npos) {
            in_ports = true;
        } else if (in_ports) {
            if (line == ");") {
                in_ports = false;
            } else if (!line.empty()) {
                port_lines.push_back(i);
            }
        }
    }

    for (size_t idx = 0; idx < port_lines.size(); ++idx) {
        auto& line = lines[port_lines[idx]];
        auto code = code_part(line);
        if (code.size() > 0 && code.back() == ',') {
            if (idx + 1 == port_lines.size()) {
                auto const pos = line.find(',');
                line.erase(pos, 1);
            }
        } else if (idx + 1 < port_lines.size()) {
            auto const comment = line.find("//");
            line.insert((comment == std::string::npos) ? line.size() : comment, ",");
        }
    }

    std::set<std::string> outputs(port_names.begin(), port_names.end());
    std::regex const wire_re(R"(^(\s*)wire signed \[[^\]]*\] (\w+) = (.*);\s*$)");
    std::smatch match;
    for (auto& line : lines) {
        if (std::regex_match(line, match, wire_re) && outputs.contains(match[2].str())) {
            line = match[1].str() + "assign " + match[2].str() + " = " + match[3].str() + ";";
        }
    }

    std::string result;
    for (size_t i = 0; i < lines.size(); ++i) {
        result += lines[i];
        if (i + 1 < lines.size()) result += '\n';
    }
    return result;
}
