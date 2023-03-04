#include "Vector.h"
#include <bits/stdc++.h>

class FourVector : public Vector
{
public:
    FourVector(double x, double y, double z, double e) : Vector(x, y, z), e(e) {}

    double getE() const
    {
        return e;
    }

    double getM() const
    {                                                               // The function calculates the mass by using the energy (e) and the magnitude (getMag()) of the three-dimensional momentum vector.
                                                                    // It uses the mathematical formula for calculating the magnitude of the fourth dimension (time-like) of a four-vector 
                                                                    // in the Minkowski spacetime, where the square of the magnitude of a four-vector is equal to the negative of the square of 
        return std::sqrt(e * e - getMag() * getMag());              // its invariant mass (m^2 = - E^2 + p^2, where E is energy, p is momentum, and c=1 is used for the speed of light). 
                                                                    // Therefore, taking the square root of e^2 - |p|^2 gives the magnitude of the fourth component, which is the mass of the 
                                                                    // four-vector.
        
    }

    double getPt() const
    {
                                                                    // The function returns the transverse momentum (pT) of the object represented by the class instance. 
                                                                    // The x and y are the x and y components of the momentum vector, respectively, and std::sqrt() is a standard C++ function that 
        return std::sqrt(x * x + y * y);                            // calculates the square root of a value. The calculation of pT involves taking the square root of the sum of the squares of 
                                                                    // the x and y components of the momentum vector.
    }

    double getPx() const
    {
        return x;                                                   // This function gets the X-coordination of the particle
    }

    double getPy() const
    {
        return y;                                                   // This function gets the Y-coordination of the particle
    }

    double getPz() const
    {
        return z;                                                   // This function gets the Z-coordination of the particle
    }

    double getTheta() const
    {
        double p = getMag();                                        // The function calculates the angle theta of the vector with respect to the z-axis, which is the polar angle. 
        return std::acos(z / p);                                    // It first calculates the magnitude (p) of the vector using the getMag() function, and then it calculates theta 
                                                                    // as the inverse cosine of the z-component of the vector divided by its magnitude. The result is in radians.
    }

    double getPhi() const
    {
        double phi = std::atan2(y, x);                              // The function is returning the azimuthal angle phi of a FourVector, which is the angle between the x-axis and the 
        if (phi < 0)                                                // projection of the vector onto the x-y plane. The function calculates the angle using the std::atan2 function, 
        {                                                           // which takes both the x and y components of the vector as arguments, and returns the angle in radians. 
            phi += 2 * M_PI;                                        // The if-statement then ensures that the returned angle is between 0 and 2*pi, since std::atan2 returns angles 
        }                                                           // between -pi and pi.
        return phi;                                                 
    }

private:
    double e;
};