#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include <cmath>


class System {
public:
    System();
    System(const char* filename, unsigned int& N, double& L);

    void read_xyz_file(const char* filename, unsigned int& N, double& L);
    void calculateNumNeighboursSerial(double rc);
    void calculateNumNeighboursSerialOptimised(double rc);
    void calculateNumNeighboursSerialCellList(double rc);
    void calculateNumNeighboursMPI(double rc, int world_rank, int world_size);
    void calculateNumNeighboursMPIOptimised(double rc, int world_rank, int world_size);
    void calculateNumNeighboursMPIBalancedOptimised(double rc, int world_rank, int world_size);
    void calculateNumNeighboursMPICellList(double rc, int world_rank, int world_size);
    void calculateNumNeighboursOpenMPStatic(double rc);
    void calculateNumNeighboursOpenMPDynamic(double rc);
    void calculateNumNeighboursOpenMPCellList(double rc);

    std::vector< std::vector<double> > GetPosition(); 
    std::vector<int> GetNumNeighbours();
private:
    std::vector< std::vector<double> > positions;
    std::vector<int> numNeighbours;
    unsigned int numParticles;
    double length;

    std::vector< std::vector< std::vector< std::vector<int> > > > createCellList(double rc);
    void calculateNeighboursPerCell(std::vector <int> currentCell, int particleIndex, double rc);
    void calculateNeighboursPerCellOpenMP(std::vector <int> currentCell, int particleIndex, double rc);
};
