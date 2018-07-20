// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_OPTIONS_HPP
#define CFILES_OPTIONS_HPP

#include <string>
#include <vector>
#include <unordered_map>

#include <chemfiles/external/optional.hpp>

#include "Errors.hpp"
#include "utils.hpp"

namespace options {

/// Possible types for the options
struct Type {
public:
    enum type_t {
        String  = 1 << 0,
        Bool    = 1 << 1,
        Vector  = 1 << 2,
        Double  = 1 << 3,
    };

    Type(type_t type): type_(type) {}

    Type(const Type&) = default;
    Type(Type&&) = default;
    Type& operator=(const Type&) = default;
    Type& operator=(Type&&) = default;

    /// Implicit convertion to the underlying enum type
    operator type_t() const {
        return type_;
    }

    /// String representation of the type
    std::string str() const {
        switch (type_) {
        case Bool:
            return "bool";
        case Double:
            return "double";
        case String:
            return "string";
        case Vector:
            return "vector";
        }
    }

private:
    type_t type_;
};


/// Possible values for an option, represented as a tagged union
class Value {
public:
    using bool_t = bool;
    using double_t = double;
    using string_t = std::string;
    using vector_t = std::vector<double>;

    Value(bool value): bool_(std::move(value)), type_(Type::Bool) {}
    Value(double value): double_(std::move(value)), type_(Type::Double) {}
    Value(const char* value): string_(value), type_(Type::String) {}
    Value(char* value): string_(value), type_(Type::String) {}
    Value(string_t value): string_(std::move(value)), type_(Type::String) {}
    Value(vector_t value): vector_(std::move(value)), type_(Type::Vector) {}

    ~Value() {
        switch (type_) {
        case Type::Bool:
        case Type::Double:
            // nothing to do
            break;
        case Type::String:
            string_.~string_t();
            break;
        case Type::Vector:
            vector_.~vector_t();
            break;
        }
    }

    Value(const Value& other): Value(false) {
        *this = other;
    }

    Value& operator=(const Value& other) {
        this->~Value();
        type_ = other.type_;
        switch (type_) {
        case Type::Bool:
            new (&this->bool_) bool_t(other.bool_);
            break;
        case Type::Double:
            new (&this->double_) double_t(other.double_);
            break;
        case Type::String:
            new (&this->string_) string_t(other.string_);
            break;
        case Type::Vector:
            new (&this->vector_) vector_t(other.vector_);
            break;
        }
        return *this;
    }

    Value(Value&& other): Value(false) {
        *this = std::move(other);
    }

    Value& operator=(Value&& other) {
        this->~Value();
        type_ = other.type_;
        switch (type_) {
        case Type::Bool:
            new (&this->bool_) bool_t(std::move(other.bool_));
            break;
        case Type::Double:
            new (&this->double_) double_t(std::move(other.double_));
            break;
        case Type::String:
            new (&this->string_) string_t(std::move(other.string_));
            break;
        case Type::Vector:
            new (&this->vector_) vector_t(std::move(other.vector_));
            break;
        }
        return *this;
    }

    /// Get the boolean in this `Value`, if the type is `Type::Bool`. Else
    /// throws an `OptionError`.
    bool_t as_bool() const {
        if (type_ == Type::Bool) {
            return bool_;
        } else {
            throw OptionError("Can not get a bool value out of a " + type_.str() + " value");
        }
    }

    /// Get the double in this `Value`, if the type is `Type::Double`. Else
    /// throws an `OptionError`.
    double_t as_double() const {
        if (type_ == Type::Double) {
            return double_;
        } else {
            throw OptionError("Can not get a double value out of a " + type_.str() + " value");
        }
    }

    /// Get the string in this `Value`, if the type is `Type::String`. Else
    /// throws an `OptionError`.
    const string_t& as_string() const {
        if (type_ == Type::String) {
            return string_;
        } else {
            throw OptionError("Can not get a string value out of a " + type_.str() + " value");
        }
    }

    /// Get the vector in this `Value`, if the type is `Type::Vector`. Else
    /// throws an `OptionError`.
    const vector_t& as_vector() const {
        if (type_ == Type::Vector) {
            return vector_;
        } else {
            throw OptionError("Can not get a vector value out of a " + type_.str() + " value");
        }
    }

    operator bool_t() const {
        return this->as_bool();
    }

    operator double_t() const {
        return this->as_double();
    }

    operator const string_t&() const {
        return this->as_string();
    }

    operator const vector_t&() const {
        return this->as_vector();
    }

    /// Get the type of this `Value`
    Type type() const {
        return type_;
    }

private:
    // Data stprage as an union
    union {
        bool bool_;
        double double_;
        string_t string_;
        vector_t vector_;
    };
    // Type tag
    Type type_;
};


/// A single option to be passed to a programm
class Option {
public:
    /// Create a new option with the given long `name`
    Option(std::string name): name_(std::move(name)) {
        if (name_.size() < 1) {
            throw OptionError("long option is too short");
        }
    }

    /// Set the default value of this option to `value`.
    ///
    /// `value` should have the same type as the expected type of this option
    Option&& set_default(Value value) && {
        if (value.type() & type_) {
            default_ = std::move(value);
        } else {
            throw OptionError("Invalid default of type " + value.type().str() +" for option " + name_);
        }
        return std::move(*this);
    }

    /// Get the default value associated with this option, if any.
    chemfiles::optional<Value> get_default() const {
        return default_;
    }

    /// Get the name of this option
    const std::string& name() const {
        return name_;
    }

    Option&& short_name(char short_name) && {
        short_name_ = short_name;
        return std::move(*this);
    }

    /// Get the optional short name of this option, or `'\0'`
    char short_name() const {
        return short_name_;
    }

    Option&& argument(std::string argument) && {
        argument_ = std::move(argument);
        return std::move(*this);
    }

    /// Get the name of the argument for this option
    const std::string& argument() const {
        return argument_;
    }

    Option&& description(std::string description) && {
        description_ = std::move(description);
        return std::move(*this);
    }

    /// Get the description of this option
    const std::string& description() const {
        return description_;
    }

    Option&& type(int type) && {
        type_ = type;
        return std::move(*this);
    }

    int type() const {
        return type_;
    }

private:
    char short_name_ = '\0';
    std::string name_;
    std::string argument_;
    std::string description_;
    int type_ = 0;
    chemfiles::optional<Value> default_ = chemfiles::nullopt;
};

} // namespace options


/// A list of options associated with a command
class Options {
public:
    template<typename ... Opts>
    Options(std::string brief, std::string author, std::string description, Opts && ... options):
        brief_(std::move(brief)),
        author_(std::move(author)),
        description_(std::move(description)),
        available_options_({std::forward<Opts>(options)...})
    {
        clean_description();
        validate_options();
    }

    const std::string& brief() const {
        return brief_;
    }

    const std::string& author() const {
        return author_;
    }

    const std::string& description() const {
        return description_;
    }

    /// Set the option with the given `name` to the specified `value`
    void set(const std::string& name, options::Value value);

    /// Parse command line arguments to populate the expected options
    void parse(int argc, const char* argv[]);

    /// Get the full help message with long description and arguments
    std::string help() const;

    /// Get the full list of available options for this command
    const std::vector<options::Option>& options() const {
        return available_options_;
    }

    /// Add a positional arguments to this command
    void add_positional(std::string positional) {
        positionals_.emplace_back(std::move(positional));
    }

    /// Get the parsed list of positional arguments for this command
    const std::vector<std::string>& positionals() const {
        return positionals_;
    }

    /// Get the option with the given `name`, set in `parse` or `set`.
    chemfiles::optional<const options::Value&> operator[](const std::string& name);

private:
    std::string arguments_help() const;
    void clean_description();
    void validate_options() const;

    /// Single line brief description
    std::string brief_;
    /// Author of the command
    std::string author_;
    /// Multi-line long description
    std::string description_;
    /// Available options for the command
    std::vector<options::Option> available_options_;

    /// Actual options, set in parse or set.
    std::unordered_map<std::string, options::Value> options_;

    /// Actual positional arguments
    std::vector<std::string> positionals_;
};


#endif
