#include <catch.hpp>
#include <chemfiles.hpp>

#include "Axis.hpp"
#include "Errors.hpp"

static bool roughly(const Vector3D& a, const Vector3D& b, double eps=1e-12) {
    return (fabs(a[0] - b[0]) < eps) && (fabs(a[1] - b[1]) < eps) && (fabs(a[2] - b[2]) < eps);
}

TEST_CASE("Axis") {
    SECTION("Constructor") {
        auto axis = Axis(1, 0, 0, Axis::Linear);
        CHECK(axis.vector() == Vector3D(1, 0, 0));
        CHECK(axis.is_linear());

        axis = Axis(0, 1, 0, Axis::Radial);
        CHECK(axis.vector() == Vector3D(0, 1, 0));
        CHECK(axis.is_radial());

        axis = Axis(2, 0, 0, Axis::Linear);
        CHECK(axis.vector() == Vector3D(1, 0, 0));

        axis = Axis(-1, 1, 2, Axis::Linear);
        CHECK(axis.vector() == Vector3D(-1/sqrt(6), 1/sqrt(6), 2/sqrt(6)));

        CHECK_THROWS_AS(Axis(0, 0, 0, Axis::Linear), CFilesError);
    }

    SECTION("Named axis") {
        auto axis = Axis::parse("x", Axis::Linear);
        CHECK(axis.vector() == Vector3D(1, 0, 0));

        axis = Axis::parse("X", Axis::Linear);
        CHECK(axis.vector() == Vector3D(1, 0, 0));

        axis = Axis::parse("y", Axis::Linear);
        CHECK(axis.vector() == Vector3D(0, 1, 0));

        axis = Axis::parse("Y", Axis::Linear);
        CHECK(axis.vector() == Vector3D(0, 1, 0));

        axis = Axis::parse("z", Axis::Linear);
        CHECK(axis.vector() == Vector3D(0, 0, 1));

        axis = Axis::parse("Z", Axis::Linear);
        CHECK(axis.vector() == Vector3D(0, 0, 1));
    }

    SECTION("Parse") {
        auto axis = Axis::parse("1:1:1", Axis::Linear);
        CHECK(roughly(axis.vector(), Vector3D(1/sqrt(3), 1/sqrt(3), 1/sqrt(3))));

        axis = Axis::parse("1:1:2", Axis::Linear);
        CHECK(roughly(axis.vector(), Vector3D(1/sqrt(6), 1/sqrt(6), 2/sqrt(6))));

        axis = Axis::parse("-1:1:2", Axis::Linear);
        CHECK(roughly(axis.vector(), Vector3D(-1/sqrt(6), 1/sqrt(6), 2/sqrt(6))));
    }

    SECTION("Parse errors") {
        auto BAD_AXIS = {
            "xy",
            "top",
            "0,1,2",
            "0:1,2",
            "0:1:2:3",
        };

        for (auto axis: BAD_AXIS) {
            CHECK_THROWS_AS(Axis::parse(axis, Axis::Linear), CFilesError);
        }
    }

    SECTION("Axis to string") {
        CHECK(Axis(1, 0, 0, Axis::Linear).str() == "x");
        CHECK(Axis(0, 1, 0, Axis::Linear).str() == "y");
        CHECK(Axis(0, 0, 1, Axis::Linear).str() == "z");
        CHECK(Axis(1, 4, -3, Axis::Linear).str() == "(0.196116, 0.784465, -0.588348)");
    }

    SECTION("Projections") {
        auto x = Axis(1, 0, 0, Axis::Linear);
        CHECK(x.projection(Vector3D(4, 0, 0)) == 4);
        CHECK(x.projection(Vector3D(-1.3, 0, 0)) == -1.3);
        CHECK(x.projection(Vector3D(3, 456, 28)) == 3);

        auto z = Axis(0, 0, 1, Axis::Radial);
        CHECK(z.projection(Vector3D(4, 0, 0)) == 4);
        CHECK(z.projection(Vector3D(-1.3, 0, 0)) == 1.3);
        CHECK(z.projection(Vector3D(3, 4, 28)) == sqrt(3 * 3 + 4 * 4));
    }
}
