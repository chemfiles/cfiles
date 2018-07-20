// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "Options.hpp"
#include "utils.hpp"

#include <fmt/format.h>

#include <unordered_set>
#include <cassert>

using namespace options;

template <typename... Args>
static OptionError option_error(const char *format, const Args & ... arguments) {
    return OptionError(fmt::format(format, arguments...));
}

/// Wrap a given string at 80 chars, ignoring single new line and repeated
/// spaces in the input string. The output is indented by `shift` spaces.
///
/// This function assumes input is in the ASCII block
class Wrapper {
public:
    Wrapper(const std::string& input, unsigned shift): input_(input), shift_(shift) {}

    std::string wrap() {
        std::string output = std::string(shift_, ' ');
        for (char c: input_) {
            if (line_.size() >= 78) {
                emit_line(output);
            }

            if (c == '\n' || c == ' ' || c == '\t') {
                if (!previous_was_space_) {
                    emit_word(output);
                }
                previous_was_space_ = true;
                continue;
            } else {
                previous_was_space_ = false;
                word_ += c;
                continue;
            }
        }
        emit_word(output);
        output += line_;
        line_.clear();

        return output;
    }

private:
    void emit_line(std::string& output) {
        output += line_ + "\n";
        output += std::string(shift_, ' ');
        line_.clear();
    }

    void emit_word(std::string& output) {
        if (line_.size() + word_.size() >= 78) {
            emit_line(output);
        }
        line_ += word_ + " ";
        word_.clear();
    }

    const std::string& input_;
    unsigned shift_;
    std::string line_;
    std::string word_;
    bool previous_was_space_;
};


/// Parser using command line arguments (argv/argc) to fill an Options instance
/// with the right values.
class OptionsParser {
public:
    OptionsParser(int argc, const char* argv[]): arguments_(argv, argv + argc) {}

    void parse(Options& options) {
        options_ = &options;

        for (const auto& option: options_->options()) {
            auto default_value = option.get_default();
            if (default_value) {
                options_->set(option.name(), std::move(*default_value));
            }
        }

        while (!done()) {
            if (peek() == "--") {
                break;  // Remaining arguments are positional arguments
            }

            if (peek()[0] == '-') {
                if (peek().size() >= 2 && peek()[1] == '-') {
                    // long option
                    auto name = advance().substr(2);
                    auto equal_pos = name.find('=');
                    if (equal_pos != std::string::npos) {
                        // --foo=bar
                        auto& option = get_option(name.substr(0, equal_pos));
                        parse_option(option, name.substr(equal_pos + 1));
                    } else {
                        // --foo or --foo bar
                        auto& option = get_option(name);
                        if (option.type() & Type::Bool) {
                            options_->set(option.name(), true);
                        } else {
                            parse_option(option, advance());
                        }
                    }
                } else {
                    // short option
                    auto name = advance();
                    auto& option = get_short_option(name[1]);
                    if (name.size() > 2) {
                        // -fbar
                        parse_option(option, name.substr(2));
                    } else {
                        // -f or -f bar
                        if (option.type() & Type::Bool) {
                            options_->set(option.name(), true);
                        } else {
                            parse_option(option, advance());
                        }
                    }
                }
            } else {
                positional();
            }
        }

        // Any remaining argument is a positional argument
        while (!done()) {
            positional();
        }
    }

private:
    void positional() {
        options_->add_positional(advance());
    }

    void parse_option(const Option& option, const std::string& value) {
        if (value == "") {
            throw option_error("expected a value for '--{}'", option.name());
        }
        if (option.type() & Type::Double) {
            try {
                options_->set(option.name(), string2double(value));
                return;
            } catch (const CFilesError&) {
                // Do nothing and continue in case we can parse this as a string
            }
        } else if (option.type() & Type::Vector) {
            try {
                auto splitted = split(value, ':');
                auto vector = std::vector<double>();
                for (auto& v: splitted) {
                    vector.push_back(string2double(v));
                }
                options_->set(option.name(), vector);
                return;
            } catch (const CFilesError&) {
                // Do nothing and continue in case we can parse this as a string
            }
        } else if (option.type() & Type::String) {
            options_->set(option.name(), value);
            return;
        }
        throw option_error("'{}' does not make sense for '--{}'", value, option.name());
    }


    const Option& get_option(const std::string& name) {
        for (auto& option: options_->options()) {
            if (option.name() == name) {
                return option;
            }
        }
        throw option_error("unknown option --{}", name);
    }

    const Option& get_short_option(char name) {
        for (auto& option: options_->options()) {
            if (option.short_name() == name) {
                return option;
            }
        }
        throw option_error("unknown option -{}", name);
    }


    bool done() const {
        return current_ >= arguments_.size();
    }

    const std::string& advance() {
        if (done()) {
            static std::string EMPTY;
            return EMPTY;
        } else {
            current_++;
            return arguments_[current_ - 1];
        }
    }

    const std::string& peek() const {
        if (done()) {
            static std::string EMPTY = "";
            return EMPTY;
        } else {
            return arguments_[current_];
        }
    }

    std::vector<std::string> arguments_;
    Options* options_ = nullptr;
    size_t current_ = 0;
};

void Options::clean_description() {
    // Find the first newline and the following number of spaces
    auto it = description_.find('\n');
    if (it == std::string::npos) {
        return;
    }

    size_t nspaces = 0;
    for (auto i = it + 1; i < description_.size(); i++) {
        if (description_[i] == ' ') {
            nspaces += 1;
        } else {
            break;
        }
    }

    // Remove the global indentation
    auto indent = "\n" + std::string(nspaces, ' ');
    it = 0;
    while (true) {
        it = description_.find(indent, it);
        if (it != std::string::npos) {
            description_.replace(it, indent.size(), "\n");
            it -= indent.size();
        } else {
            break;
        }
    }
}

void Options::validate_options() const {
    std::unordered_set<std::string> names;
    std::unordered_set<char> short_names;
    for (auto& option: available_options_) {
        auto it = names.insert(option.name());
        if (!it.second) {
            throw option_error("The option {} was specified twice", option.name());
        }

        if (option.short_name() != '\0') {
            auto it = short_names.insert(option.short_name());
            if (!it.second) {
                throw option_error("The short option name {} was specified twice", option.short_name());
            }
        }

        if (option.type() == 0) {
            throw option_error("Please specify a type for the '--{}' option", option.name());
        }

        if (option.description() == "") {
            throw option_error("Please specify a description for the '--{}' option", option.name());
        }
    }
}

void Options::set(const std::string& name, Value value) {
    for (auto& option: available_options_) {
        if (option.name() == name) {
            if (option.type() & value.type()) {
                auto it = options_.emplace(name, value);
                if (!it.second) {
                    it.first->second = value;
                }
                return;
            } else {
                throw option_error("Invalid type {} for '--{}' option", value.type().str(), option.name());
            }
        }
    }
    throw option_error("Could not find an option called {}", name);
}

chemfiles::optional<const Value&> Options::operator[](const std::string& name) {
    auto it = options_.find(name);
    if (it != options_.end()) {
        return it->second;
    }
    return chemfiles::nullopt;
}

void Options::parse(int argc, const char* argv[]) {
    auto parser = OptionsParser(argc, argv);
    parser.parse(*this);
}

std::string Options::help() const {
    std::string output;
    output += brief_ + "\n";
    output += full_version() + "\n";
    output += author_ + "\n\n";
    output += description_ + "\n\n";
    output += arguments_help() + "\n";
    return output;
}

std::string Options::arguments_help() const {
    fmt::memory_buffer output;
    fmt::format_to(output, "Options:\n");
    for (const auto& option: available_options_) {
        fmt::format_to(output, "  ");
        if (option.short_name() != '\0') {
            fmt::format_to(output, "-{}", option.short_name());
            if (option.argument() != "") {
                fmt::format_to(output, " <{}>", option.argument());
            }
            fmt::format_to(output, ", ");
        }
        fmt::format_to(output, "--{}", option.name());
        if (option.argument() != "") {
            fmt::format_to(output, "=<{}>", option.argument());
        };
        auto description = option.description();

        auto default_value = option.get_default();
        if (default_value && default_value->type() == Type::String) {
            description += fmt::format(" [default: {}]", default_value->as_string());
        } else if (default_value && default_value->type() == Type::Double) {
            description += fmt::format(" [default: {}]", default_value->as_double());
        }
        fmt::format_to(output, "\n{}\n\n", Wrapper(description, 6).wrap());
    }
    return fmt::to_string(output).substr(0, output.size() - 2);
}
