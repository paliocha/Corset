// Lightweight progress / formatting utilities for console output.
// Header-only — no .cc file needed.  All functions are inline.
//
// Features:
//   - TTY detection (cached once, thread-safe)
//   - ANSI color wrappers (no-op when stdout is not a TTY)
//   - Human-friendly formatting for durations, counts, and rates
//   - ProgressLine: overwriting single-line progress for serial stages
//   - Stage banners for visual structure
//
// Dependencies: <string>, <cstdio>, <cstdint>, <cmath>, <omp.h>, <unistd.h>
//
// Author: Martin Paliocha, martin.paliocha@nmbu.no
// Created 23 February 2026

#pragma once

#include <string>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <omp.h>
#include <unistd.h>

namespace progress {


// ── TTY detection ────────────────────────────────────────────────────

// Cached on first call.  C++11 guarantees thread-safe initialization
// of function-local statics — safe to call from any thread.
inline bool is_tty() {
    static const bool tty = isatty(fileno(stdout));
    return tty;
}


// ── ANSI color wrappers ─────────────────────────────────────────────
// Each returns the plain string when stdout is not a TTY.

namespace ansi {

inline std::string bold(const std::string &s) {
    return is_tty() ? "\033[1m" + s + "\033[0m" : s;
}
inline std::string dim(const std::string &s) {
    return is_tty() ? "\033[2m" + s + "\033[0m" : s;
}
inline std::string green(const std::string &s) {
    return is_tty() ? "\033[32m" + s + "\033[0m" : s;
}
inline std::string yellow(const std::string &s) {
    return is_tty() ? "\033[33m" + s + "\033[0m" : s;
}
inline std::string cyan(const std::string &s) {
    return is_tty() ? "\033[36m" + s + "\033[0m" : s;
}
inline std::string red(const std::string &s) {
    return is_tty() ? "\033[31m" + s + "\033[0m" : s;
}

} // namespace ansi


// ── Duration formatting ─────────────────────────────────────────────

inline std::string format_duration(double secs) {
    if (secs < 0.01) return "0.0s";
    if (secs < 60.0) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.1fs", secs);
        return buf;
    }
    if (secs < 3600.0) {
        int m = static_cast<int>(secs) / 60;
        int s = static_cast<int>(secs) % 60;
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%dm%02ds", m, s);
        return buf;
    }
    int h = static_cast<int>(secs) / 3600;
    int m = (static_cast<int>(secs) % 3600) / 60;
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%dh%02dm", h, m);
    return buf;
}

inline std::string format_eta(double remaining_secs) {
    if (remaining_secs < 0) return "~?";
    return "~" + format_duration(remaining_secs);
}


// ── Count / rate formatting ─────────────────────────────────────────

inline std::string format_count(int64_t n) {
    if (n < 1000) return std::to_string(n);
    if (n < 1000000) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.1fK", n / 1e3);
        return buf;
    }
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%.2fM", n / 1e6);
    return buf;
}

inline std::string format_rate(double per_sec) {
    if (per_sec < 1.0) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.2f/s", per_sec);
        return buf;
    }
    return format_count(static_cast<int64_t>(per_sec)) + "/s";
}


// ── Progress bar (UTF-8 block elements) ─────────────────────────────

inline std::string progress_bar(double fraction, int width = 20) {
    int filled = static_cast<int>(fraction * width);
    if (filled > width) filled = width;
    if (filled < 0) filled = 0;
    std::string bar;
    for (int i = 0; i < filled; i++)     bar += "\xe2\x96\x88";  // U+2588 FULL BLOCK
    for (int i = filled; i < width; i++) bar += "\xe2\x96\x91";  // U+2591 LIGHT SHADE
    return bar;
}


// ── ProgressLine: overwriting progress for serial stages ────────────
//
// Two modes:
//   Known-total:   ProgressLine("label", total)  — bar + % + rate + ETA
//   Unknown-total: ProgressLine("label")          — count + rate
//
// TTY:     \r-overwriting updates (throttled to 0.25s)
// Non-TTY: line-based updates (throttled to 5s)
// finish() always prints a final summary with \n.

class ProgressLine {
    std::string label_;
    int64_t total_;           // -1 = unknown total (count-up mode)
    int64_t current_ = 0;
    double start_;
    double last_print_ = 0;
    bool tty_;

public:
    ProgressLine(const std::string &label, int64_t total)
        : label_(label), total_(total), start_(omp_get_wtime()), tty_(is_tty()) {}

    explicit ProgressLine(const std::string &label)
        : label_(label), total_(-1), start_(omp_get_wtime()), tty_(is_tty()) {}

    void update(int64_t current) {
        current_ = current;
        double now = omp_get_wtime();
        double interval = tty_ ? 0.25 : 5.0;
        if (now - last_print_ < interval) return;
        last_print_ = now;
        print(false);
    }

    void finish() {
        print(true);
    }

private:
    void print(bool final) {
        double elapsed = omp_get_wtime() - start_;
        double rate = (current_ > 0 && elapsed > 0.001)
                          ? static_cast<double>(current_) / elapsed : 0;

        if (tty_) {
            std::string line;
            if (total_ > 0) {
                double frac = static_cast<double>(current_) / total_;
                double eta = (rate > 0) ? (total_ - current_) / rate : -1;
                char pct[8];
                std::snprintf(pct, sizeof(pct), "%5.1f%%", frac * 100);
                line = "[" + label_ + "] " + progress_bar(frac)
                     + " " + pct + "  " + format_count(current_)
                     + "  " + format_rate(rate);
                if (!final && eta > 0) line += "  " + format_eta(eta);
            } else {
                line = "[" + label_ + "] " + format_count(current_)
                     + "  " + format_rate(rate);
            }
            if (final) line += "  " + format_duration(elapsed);

            // Pad to overwrite previous line remnants
            while (line.size() < 80) line += ' ';

            if (final) std::printf("\r%s\n", line.c_str());
            else       std::printf("\r%s", line.c_str());
            std::fflush(stdout);
        } else {
            // Non-TTY: single line on finish, periodic lines otherwise
            if (final) {
                std::printf("[%s] %s in %s (%s)\n",
                            label_.c_str(), format_count(current_).c_str(),
                            format_duration(elapsed).c_str(),
                            format_rate(rate).c_str());
            } else {
                std::printf("[%s] %s (%s)\n",
                            label_.c_str(), format_count(current_).c_str(),
                            format_rate(rate).c_str());
            }
            std::fflush(stdout);
        }
    }
};


// ── Stage banners ───────────────────────────────────────────────────

inline void print_banner(const std::string &text) {
    // "── text ──────────..."  (UTF-8 box-drawing horizontal line)
    std::string line = "\xe2\x94\x80\xe2\x94\x80 " + text + " ";
    while (line.size() < 60) line += "\xe2\x94\x80";
    std::printf("\n%s\n", ansi::bold(line).c_str());
    std::fflush(stdout);
}

inline void print_separator() {
    std::string sep;
    for (int i = 0; i < 57; i++) sep += "\xe2\x94\x80";
    std::printf("%s\n", ansi::dim(sep).c_str());
    std::fflush(stdout);
}


} // namespace progress
