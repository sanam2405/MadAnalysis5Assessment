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