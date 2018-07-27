// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <DefaultOptions.hpp>

Option OPTION_OUTPUT(const std::string& extension) {
    return Option("output")
        .short_name('o')
        .argument("file")
        .description(
            R"(write result to <file>. This default to the input trajectory file
            name with the )" + extension + " extension."
        );
}

Option OPTION_FORMAT = Option("format")
    .argument("format")
    .description("force the input file format to be <format>");

Option OPTION_TOPOLOGY = Option("topology")
    .short_name('t')
    .argument("path")
    .description("alternative topology file for the input");

Option OPTION_TOPOLOGY_FORMAT = Option("topology-format")
    .argument("format")
    .description("use <format> as format for the topology file");

Option OPTION_GUESS_BONDS = Option("guess-bonds")
    .description("guess the bonds in the input");

Option OPTION_CELL = Option("cell")
    .short_name('c')
    .argument("cell")
    .description(
        R"(alternative unit cell. <cell> format is one of <a:b:c:α:β:γ> or
        <a:b:c> or <a>. 'a', 'b' and 'c' are in angstroms; 'α', 'β', and 'γ' are
        in degrees.)"
    );

Option OPTION_STEPS = Option("steps")
    .argument("steps")
    .description(
        R"(steps to use from the input. <steps> format is
        <start>:<end>[:<stride>] with <start>, <end> and <stride> optional. The
        used steps goes from <start> to <end> (excluded) by steps of <stride>.
        The default values are 0 for <start>, the number of steps for <end> and
        1 for <stride>)"
    );
