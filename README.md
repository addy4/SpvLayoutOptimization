# SPV_Layout_Optimization
* A genetic algorithm has been designed for optimizing the cost of solar PV power plants. This algorithm generates an optimal SPV layout for a given SPV plant. 
The input to this algorithm includes data about the number of solar PV arrays, number of battery-banks to be placed and the locations of the SPV arrays in this plant.

* This repository consists of input test instances designed for the aforesaid genetic algorithm. Also, the generation code for designing these test instances has been
shared
* Following commands need to be invoked for creating the test instance:
  - cd generateTest
  - vim testcases.cpp (and make changes in X, Y, N, P values as per requirement) 
  - g++ -std=c++17 testcases.cpp -o TestC 
  - ./TestC > ../testInstances/CaseP_N.rt 
 
