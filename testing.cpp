#include <catch2/catch_all.hpp> // Include the Catch2 library
#include <catch2/catch_approx.hpp>
#include "FourVector.h" // Include the header file for the FourVector class

TEST_CASE("FourVector basic properties", "[FourVector]")
{
    // Test the constructor and getters
    FourVector v(1, 2, 3, 4);
    REQUIRE(v.getX() == 1);
    REQUIRE(v.getY() == 2);
    REQUIRE(v.getZ() == 3);
    REQUIRE(v.getE() == 4);

    // Test the mass calculation
    REQUIRE(v.getM() == Approx(2.236068));

    // Test the transverse momentum calculation
    REQUIRE(v.getPt() == Approx(2.236068));

    // Test the momentum component getters
    REQUIRE(v.getPx() == 1);
    REQUIRE(v.getPy() == 2);
    REQUIRE(v.getPz() == 3);

    // Test the theta calculation
    REQUIRE(v.getTheta() == Approx(0.982793));

    // Test the phi calculation
    REQUIRE(v.getPhi() == Approx(1.107149));
}
