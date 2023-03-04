
#include <cmath>

class Vector
{
public:
    // Constructor that takes in x, y, and z components of the vector
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    // Getter functions to access x, y, and z components of the vector
    double getX() const
    {
        return x;
    }

    double getY() const
    {
        return y;
    }

    double getZ() const
    {
        return z;
    }

    // Calculates and returns the magnitude of the vector
    double getMag() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Calculates and returns a unit vector in the same direction as the original vector
    Vector getUnit() const
    {
        double mag = getMag();
        return Vector(x / mag, y / mag, z / mag);
    }

    // Calculates and returns the dot product of the vector with another vector
    double dot(const Vector &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    // Calculates and returns the cross product of the vector with another vector
    Vector cross(const Vector &other) const
    {
        double x_cross = y * other.z - z * other.y;
        double y_cross = z * other.x - x * other.z;
        double z_cross = x * other.y - y * other.x;
        return Vector(x_cross, y_cross, z_cross);
    }

protected:
    // Private member variables to store the x, y, and z components of the vector
    double x, y, z;
};