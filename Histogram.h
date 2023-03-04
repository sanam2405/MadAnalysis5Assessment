#include <bits/stdc++.h>

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