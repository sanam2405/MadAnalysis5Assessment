#include <iostream>     // for basic input and output
#include <vector>       // for using vector data structure
#include <cmath>        // for mathematical computations
#include <algorithm>    // for std::upper_bound
#include <numeric>      // for std::accumulate
#include <fstream>      // for file output
#include <iomanip>      // for output formatting


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

class Histogram
{
public:
    std::vector<double> data; // Vector to hold the histogram data
    std::vector<double> bins; // Vector to hold the bin edges

    // Constructor to initialize the histogram
    Histogram(int nBins, double binMin, double binMax)
    {
        bins.resize(nBins + 1); // Resize the bin vector to the desired number of bins
        for (int i = 0; i <= nBins; ++i)
        {
            bins[i] = binMin + (binMax - binMin) / nBins * i; // Calculate the bin edges
        }
        data.resize(nBins); // Resize the data vector to the desired number of bins
    }

    // Function to fill the histogram with a given value
    void fill(double value)
    {
        int bin = std::upper_bound(bins.begin(), bins.end(), value) - bins.begin() - 1; // Find the bin index where the value belongs
        if (bin >= 0 && bin < data.size()) // Check if the bin index is within the range of the data vector
        {
            ++data[bin]; // Increment the count for the corresponding bin
        }
    }

    // Function to normalize the histogram data
    void normalize()
    {
        double zero = (double) 0.0;
        double sum = std::accumulate(data.begin(), data.end(), zero); // Calculate the sum of all the histogram counts
        if (sum > 0) // Check if the sum is greater than zero to avoid division by zero
        {
            for (auto &d : data)
            {
                d /= sum; // Divide each histogram count by the total sum to obtain normalized counts
            }
        }
    }

    // Function to save the histogram data to a file
    void save(std::string filename)
    {
        std::ofstream fout(filename); // Open a file stream for writing
        fout << "BINS"
             << "\t"
             << "DATA"
             << "\n"; // Write the header line for the data file
        for (int i = 0; i < bins.size() - 1; ++i) // Loop over all the bins and write the bin center and corresponding count to the file
        {
            fout << std::setprecision(10) << 0.5 * (bins[i] + bins[i + 1]) << "\t" << std::setprecision(10) << data[i] << "\n";
        }
        fout.close(); // Close the file stream
    }
};

class Reader
{
public:
    // Constructor that takes a filename as input
    Reader(const std::string &filename)
    {
        // Open the file with the given filename
        std::ifstream file(filename);
        
        // Declare variables for x, y, z, and e
        double x, y, z, e;
        
        // Read in the data from the file
        while (file >> x >> y >> z >> e)
        {
            // Create a FourVector object from the data
            FourVector vec(x, y, z, e);
            
            // Add the FourVector object to the vector of FourVectors
            vectors.push_back(vec);
        }
        file.close();
    }
    
    // Method to return the vector of FourVectors
    std::vector<FourVector> getVectors() const
    {
        return vectors;
    }
    
private:
    std::vector<FourVector> vectors;  // Vector to store FourVectors
};


int main(int argc, char *argv[])
{
    // parsing command-line arguments
    if (argc != 6)
    {
        std::cout << "Usage: " << argv[0] << " datafile distribution_name n_bins x_min x_max\n";
        return 0;
    }
    std::string datafile = argv[1];
    std::string distname = argv[2];
    int nBins = std::atoi(argv[3]);
    double xMin = std::atof(argv[4]);
    double xMax = std::atof(argv[5]);

    // reading data file
    std::vector<FourVector> particles;
    // std::ifstream fin(datafile);
    // double x, y, z, E;
    // while (fin >> x >> y >> z >> E)
    // {
    //     particles.emplace_back(x, y, z, E);
    // }
    // fin.close();

    Reader read(datafile);
    particles = read.getVectors();

    // creating histogram
    Histogram hist(nBins, xMin, xMax);

    // filling histogram based on distribution name
    if (distname == "pT")
    {
        for (const auto &p : particles)
        {
            hist.fill(p.getPt());
        }
    }
    else if (distname == "pX")
    {
        for (const auto &p : particles)
        {
            hist.fill(p.getPx());
        }
    }
    else if (distname == "pY")
    {
        for (const auto &p : particles)
        {
            hist.fill(p.getPy());
        }
    }
    else if (distname == "pZ")
    {
        for (const auto &p : particles)
        {
            hist.fill(p.getPz());
        }
    }
    else if (distname == "energy")
    {
        for (const auto &p : particles)
        {
            hist.fill(p.getE());
        }
    }
    else if (distname == "mass")
    {
        for (const auto &p : particles)
        {
            hist.fill(p.getM());
        }
    }
    else
    {
        std::cout << "Error: unknown distribution name: " << distname << "\n";
        return 1;
    }

    // normalizing histogram and saving to file
    hist.normalize();
    std::string outfilename = distname + "_histogram.dat";
    hist.save(outfilename);

    return 0;
}
