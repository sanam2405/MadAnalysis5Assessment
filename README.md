# MadAnalysis5Assessment

## My Understanding of the Problem Statement

The problem requires to write a C++ code that reads data from a file containing four-vector data. Each four-vector consists of three momentum components and an energy component, all given in units of MeV. The code should have a reader class that reads the data and creates an appropriate data class. The data class should consist of two layers: a vector class that performs basic vector computations, and a FourVector class that inherits the vector class and adds additional capabilities such as computing mass. The code should also have an output system that generates 1D histograms for different distributions of the four-vectors, such as transverse momentum (pT), momentum in the x-direction (pX), momentum in the y-direction (pY), momentum in the z-direction (pZ), energy, and mass. The main program should take input parameters such as datafile, distribution name, number of bins, and axis limits, and generate the corresponding histogram in a two-column text file format. The code should be designed with extensive documentation and modularity, such that adding new capabilities to the FourVector class should not require significant modifications to the main program. Additionally, the code should include unit tests and be uploaded to GitHub with an appropriate workflow.

## Workflow

1. First, the required header files are included in the code, namely 'iostream', 'vector', 'cmath', 'algorithm', 'numeric', 'fstream', and 'iomanip'.

2. The Vector class is defined with the following member functions:
* A constructor that takes the x, y, and z components of a vector as arguments and initializes the corresponding private member variables.
* Three getter functions to access the private x, y, and z member variables.
* A member function getMag() that calculates and returns the magnitude of the vector using the std::sqrt() function and the Pythagorean theorem.
* A member function getUnit() that calculates and returns a unit vector in the same direction as the original vector by dividing each component of the vector by its magnitude.
* A member function dot() that calculates and returns the dot product of the vector with another vector.
* A member function cross() that calculates and returns the cross product of the vector with another vector.

3. The FourVector class is defined as a derived class of the Vector class with the following member functions:
* A constructor that takes the x, y, z, and e components of a four-vector as arguments, initializes the corresponding x, y, and z member variables using the base class constructor, and initializes the private e member variable.
* A getter function getE() to access the private e member variable.
* A member function getM() that calculates and returns the mass of the four-vector using the energy and magnitude of the three-dimensional momentum vector, and the mathematical formula for the magnitude of the fourth dimension of a four-vector in the Minkowski spacetime.
* A member function getPt() that calculates and returns the transverse momentum of the four-vector by taking the square root of the sum of the squares of the x and y components of the momentum vector.
Three getter functions getPx(), getPy(), and getPz() to access the x, y, and z components of the four-vector, respectively.
* A member function getTheta() that calculates and returns the polar angle of the vector with respect to the z-axis using the std::acos() function and the inverse cosine of the z-component of the vector divided by its magnitude.
* A member function getPhi() that calculates and returns the azimuthal angle of the four-vector using the std::atan2() function and the angle between the x-axis and the projection of the vector onto the x-y plane.

4. The Histogram class is defined with the following member variables and functions:
* Two public member variables, a std::vector<double> named data to hold the histogram data and a std::vector<double> named bins to hold the bin edges.
* A constructor that takes the minimum and maximum values of the data range and the number of bins as arguments, calculates the bin width, and initializes the bins vector.
* A member function fill() that takes a double value as an argument and adds it to the appropriate bin in the data vector using the std::upper_bound() function.
* A member function normalize() that normalizes the histogram data.
* A member function save() that takes a string filename and stores the histogram data into that file.
  
5. The Reader class is defined with the following member variables and functions
  
* A constructor that takes a file name from which the data is to be read.
* A private member varible vectors which stores the FourVectors.
* A member method which returns the vector of FourVectors.
  
  
6. The Main Method
 
* The main method is called with command-line arguments specifying the data file name, distribution name, number of bins, minimum value for x-axis, and maximum value for x-axis.
* The program checks if the number of arguments is correct (i.e., 6).
The data file name, distribution name, number of bins, minimum value for x-axis, and maximum value for x-axis are stored in their respective variables.
* An empty vector of FourVectors called "particles" is created.
* The data file is opened for reading using an input file stream.
* A loop reads in the x, y, z, and energy values from the file and creates a new FourVector object using those values. The FourVector object is then added to the "particles" vector.
* The file is closed.
* A Histogram object called "hist" is created with the specified number of bins, minimum value for x-axis, and maximum value for x-axis.
* Depending on the distribution name specified, the "hist" object is filled with values from the corresponding component of each FourVector object in the "particles" vector. If the distribution name is not recognized, an error message is displayed and the program terminates with an exit code of 1.
* The "hist" object is normalized by dividing each bin value by the total number of entries in the histogram.
* A string variable "outfilename" is created by appending "_histogram.dat" to the distribution name.
* The histogram data and bin edges are saved to a file with the name specified in "outfilename".
* The main function returns 0, indicating successful execution of the program.
