#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <chemfiles.hpp>

#include "utils.hpp"
#include "Errors.hpp"

TEST_CASE("Steps ranges") {
    auto range = steps_range::parse("10:20");
    auto result = std::vector<size_t>(range.begin(), range.end());
    auto expected = std::vector<size_t>{10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    CHECK(result == expected);

    range = steps_range::parse("10:20:2");
    result = std::vector<size_t>(range.begin(), range.end());
    expected = std::vector<size_t>{10, 12, 14, 16, 18};
    CHECK(result == expected);

    range = steps_range::parse(":20:3");
    result = std::vector<size_t>(range.begin(), range.end());
    expected = std::vector<size_t>{0, 3, 6, 9, 12, 15, 18};
    CHECK(result == expected);

    // Warning: this is an open ended range, do not create a vector from it
    range = steps_range::parse("200::5");
    auto it = range.begin();
    CHECK(*it == 200);
    ++it;
    CHECK(*it == 205);
    ++it;
    CHECK(*it == 210);
    ++it;
    CHECK(*it == 215);

    SECTION("Errors") {
        auto bad_ranges = {
            "1,2,3",
            "a:3",
            "3:a",
            "3:5:a",
            "-4:7",
            "4:-7",
            "152:3",
            "1:3:0",
        };

        for (auto& range: bad_ranges) {
            CHECK_THROWS_AS(steps_range::parse(range), CFilesError);
        }
    }
}


TEST_CASE("String to numbers") {
    CHECK(string2double("12.5") == 12.5);
    CHECK(string2long("12") == 12);

    CHECK_THROWS_AS(string2double("foo"), CFilesError);
    CHECK_THROWS_AS(string2long("foo"), CFilesError);

    CHECK_THROWS_AS(string2double("1,2"), CFilesError);
    CHECK_THROWS_AS(string2long("12,35"), CFilesError);
}

TEST_CASE("Split") {
    auto splitted = split("a,b,c,,d,  a, ", ',');
    auto expected = std::vector<std::string>{"a", "b", "c", "", "d", "  a", " "};
    CHECK(splitted == expected);
}


TEST_CASE("Parse cell") {
    auto cell = parse_cell("10");
    CHECK(cell == chemfiles::UnitCell(10));

    cell = parse_cell("10:20:30");
    CHECK(cell == chemfiles::UnitCell(10, 20, 30));

    cell = parse_cell("10:20:30:90:90:100");
    CHECK(cell == chemfiles::UnitCell(10, 20, 30, 90, 90, 100));

    SECTION("Errors") {
        auto bad_cells = {
            "1,2,3",
            "a:5:3",
            "3:a:1",
            "3:5:a",
            "4:7",
            "4:7:8:1",
            "4:7:8:1:10",
            "-7",
            "-7:4:5",
            "7:-4:5",
            "7:4:-5",
            "7:7:7:-90:90:90",
            "7:7:7:90:-90:90",
            "7:7:7:90:90:-90",
            "0",
            "0:4:5",
            "7:0:5",
            "7:4:0",
            "7:7:7:0:90:90",
            "7:7:7:90:0:90",
            "7:7:7:90:90:0",
        };

        for (auto& cell: bad_cells) {
            CHECK_THROWS_AS(parse_cell(cell), CFilesError);
        }
    }
}
