// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_UTILS_HPP
#define CFILES_UTILS_HPP

#include <string>
#include <vector>

namespace chemfiles {
    class UnitCell;
}

/// Convert a string to double
double string2double(const std::string& string);

/// Convert a string to long
long string2long(const std::string& string);

/// Get a string describing the full version of cfiles
std::string full_version();

/// Create the command description header
std::string command_header(std::string name, std::string description);

/// Split a string a delimiter
std::vector<std::string> split(const std::string& string, char delimiter);

/// Parse an unit cell string
chemfiles::UnitCell parse_cell(const std::string& string);

/// Range of steps to use from a trajectory
class steps_range {
public:
    class iterator {
       friend class steps_range;
    public:
        size_t operator*() const {return step_;}
        const iterator& operator++() {
            step_ += stride_;
            return *this;
        }
        bool operator ==(const iterator &other) const { return step_ == other.step_; }
        bool operator !=(const iterator &other) const { return step_ != other.step_; }

    protected:
        iterator(size_t start, size_t stride) : step_(start), stride_(stride) {}

    private:
        size_t step_;
        size_t stride_;
    };

    iterator begin() const {return iterator(first_, stride_);}
    iterator end() const {return iterator(last_, 0);}

    static steps_range parse(const std::string& string);
private:
    /// Starting step
    size_t first_ = 0;
    /// Last step to use
    size_t last_ = static_cast<size_t>(-1);
    /// Use one step every `stride` steps
    size_t stride_ = 1;
};

#endif
