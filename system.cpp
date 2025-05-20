#include "system.h"


System::System() {
    // default constructor
}

System::System(const char* filename, unsigned int& N, double& L) {
    // constructor that reads in the positions from a file
    read_xyz_file(filename, N, L);
}


void System::read_xyz_file(const char* filename, unsigned int& N, double& L){
    // open the file
    std::ifstream xyz_file(filename);   
    // read in the number of atoms
    xyz_file >> N;

    // read in the cell dimension
    xyz_file >> L;

    // store the number of particles and the cell dimension
    numParticles = N;
    length = L;

    // now read in the positions, ignoring the atomic species
    std::vector <double> pos = {0.0, 0.0, 0.0};
    std::string dummy; 
    for (int i=0;i<N;i++) {
        xyz_file >> dummy >> pos[0] >> pos[1] >> pos[2];
        positions.push_back(pos);           
    }

    // close the file
    xyz_file.close();
}


std::vector< std::vector<double> > System::GetPosition() {
    // return the positions
    return positions;
}

std::vector<int> System::GetNumNeighbours() {
    // return the number of neighbours
    return numNeighbours;
}


void System::calculateNumNeighboursSerial(double rc) {
    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);

    for (int i=0;i<numParticles;i++) {
        for (int j=0;j<numParticles;j++) {
            if (i!=j) {
                double r = 0.0;
                for (int k=0;k<3;k++) {
                    r += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
                }
                r = sqrt(r);
                if (r<rc) {
                    numNeighbours[i]++;
                }
            }
        }
    }
}


void System::calculateNumNeighboursSerialOptimised(double rc) {
    // square rc so that it only has to be calculated once
    double rc2 = rc*rc;

    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);

    for (int i=0;i<numParticles;i++) {
        // loop from i+1 to numParticles to reduce the number of calculations
        for (int j=i+1;j<numParticles;j++) {
            double r2 = 0.0;
            for (int k=0;k<3;k++) {
                r2 += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
            }
            if (r2<rc2) {
                numNeighbours[i]++;
                numNeighbours[j]++;
            }
        }
    }
}


std::vector< std::vector< std::vector< std::vector <int> > > > System::createCellList(double rc) {
    int numCells = (int)(length / rc);
    int cellLength = length / numCells;

    std::vector< std::vector< std::vector< std::vector<int> > > > cellList(numCells, std::vector< std::vector< std::vector<int> > >(numCells, std::vector< std::vector<int> >(numCells)));

    for (int i = 0; i < numParticles; i++) {
        int x = (int)(positions[i][0] / cellLength);
        int y = (int)(positions[i][1] / cellLength);
        int z = (int)(positions[i][2] / cellLength);

        // if the particle is at the edge of the box, it will be in the last cell
        if (x == numCells) x--;
        if (y == numCells) y--;
        if (z == numCells) z--;

        cellList[x][y][z].push_back(i);
    }

    return cellList;
}


void System::calculateNeighboursPerCell(std::vector <int> currentCell, int particleIndex, double rc2) {
    // loop over the particles in the current cell
    for (int j = 0; j < currentCell.size(); j++) {
        int neighbourIndex = currentCell[j];
        if (particleIndex != neighbourIndex) {
            double r2 = 0.0;
            for (int k = 0; k < 3; k++) {
                r2 += (positions[particleIndex][k] - positions[neighbourIndex][k]) * (positions[particleIndex][k] - positions[neighbourIndex][k]);
            }
            if (r2 < rc2) {
                numNeighbours[particleIndex]++;
            }
        }
    }
}


void System::calculateNumNeighboursSerialCellList(double rc) {
    numNeighbours.resize(numParticles, 0);

    std::vector< std::vector< std::vector< std::vector<int> > > > cellList = createCellList(rc);

    double rc2 = rc * rc;

    // loop through the 3D cell list
    for (int x = 0; x < cellList.size(); x++) {
        for (int y = 0; y < cellList[x].size(); y++) {
            for (int z = 0; z < cellList[x][y].size(); z++) {
                
                // for each cell perform the near neighbour search
                for (int i = 0; i < cellList[x][y][z].size(); i++) {
                    int currentIndex = cellList[x][y][z][i];
                    
                    // calculate the number of neighbours in all surrounding cells
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                int newX = x + dx;
                                int newY = y + dy;
                                int newZ = z + dz;

                                // calculate number of neighbours for all surrounding cells
                                if (newX >= 0 && newX < cellList.size() && newY >= 0 && newY < cellList[x].size() && newZ >= 0 && newZ < cellList[x][y].size()) {
                                    calculateNeighboursPerCell(cellList[newX][newY][newZ], currentIndex, rc2);
                                }
                            }
                        }
                    }

                }

            }
        }
    }
}

void System::calculateNumNeighboursMPI(double rc, int world_rank, int world_size) {
    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);

    // calculate the number of particles per process
    int particlesPerProcess = numParticles / world_size;
    int remainder = numParticles % world_size;

    // calculate the start and end indices for the current process
    int start = world_rank * particlesPerProcess;
    int end = (world_rank + 1) * particlesPerProcess;

    // if there is a remainder, the last process will have more particles
    if (world_rank == world_size - 1) {
        end += remainder;
    }

    for (int i=start;i<end;i++) {
        for (int j=0;j<numParticles;j++) {
            if (i!=j) {
                double r = 0.0;
                for (int k=0;k<3;k++) {
                    r += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
                }
                r = sqrt(r);
                if (r<rc) {
                    numNeighbours[i]++;
                }
            }
        }
    }

    // Use MPI_Reduce to sum all elements in the vector into the root vector
    std::vector<int> globalNumNeighbours(numParticles, 0);
    MPI_Reduce(&numNeighbours[0], &globalNumNeighbours[0], numParticles, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // if the current process is the root, set the global vector to the root vector
    if (world_rank == 0) {
        numNeighbours = globalNumNeighbours;
    }
}


void System::calculateNumNeighboursMPIOptimised(double rc, int world_rank, int world_size) {
    // square rc so that it only has to be calculated once
    double rc2 = rc*rc;

    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);

    // calculate the number of particles per process
    int particlesPerProcess = numParticles / world_size;
    int remainder = numParticles % world_size;

    // calculate the start and end indices for the current process
    int start = world_rank * particlesPerProcess;
    int end = (world_rank + 1) * particlesPerProcess;

    // if there is a remainder, the last process will have more particles
    if (world_rank == world_size - 1) {
        end += remainder;
    }

    // loop over the particles
    for (int i=start;i<end;i++) {
        // loop from i+1 to numParticles to reduce the number of calculations
        for (int j=i+1;j<numParticles;j++) {
            double r2 = 0.0;
            for (int k=0;k<3;k++) {
                r2 += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
            }
            if (r2<rc2) {
                numNeighbours[i]++;
                numNeighbours[j]++;
            }
        }
    }

    // Use MPI_Reduce to sum all elements in the vector into the root vector
    std::vector<int> globalNumNeighbours(numParticles, 0);
    MPI_Reduce(&numNeighbours[0], &globalNumNeighbours[0], numParticles, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // if the current process is the root, set the global vector to the root vector
    if (world_rank == 0) {
        numNeighbours = globalNumNeighbours;
    }
}


void System::calculateNumNeighboursMPIBalancedOptimised(double rc, int world_rank, int world_size) {
    // square rc so that it only has to be calculated once
    double rc2 = rc*rc;

    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);


    // loop over the particles in a balanced way
    for (int i = world_rank; i < numParticles; i += world_size) {
        // loop from i+1 to numParticles to reduce the number of calculations
        for (int j=i+1;j<numParticles;j++) {
            double r2 = 0.0;
            for (int k=0;k<3;k++) {
                r2 += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
            }
            if (r2<rc2) {
                numNeighbours[i]++;
                numNeighbours[j]++;
            }
        }
    }

    // Use MPI_Reduce to sum all elements in the vector into the root vector
    std::vector<int> globalNumNeighbours(numParticles, 0);
    MPI_Reduce(&numNeighbours[0], &globalNumNeighbours[0], numParticles, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // if the current process is the root, set the global vector to the root vector
    if (world_rank == 0) {
        numNeighbours = globalNumNeighbours;
    }
}


void System::calculateNumNeighboursMPICellList(double rc, int world_rank, int world_size) {
    numNeighbours.resize(numParticles, 0);
    std::vector< std::vector< std::vector< std::vector<int> > > > cellList = createCellList(rc);

    double rc2 = rc * rc;

    int counter = 0;

    // loop through the 3D cell list
    for (int x = 0; x < cellList.size(); x++) {
        for (int y = 0; y < cellList[x].size(); y++) {
            for (int z = 0; z < cellList[x][y].size(); z++) {
                
                if (counter % world_size != world_rank) {
                    counter++;
                    continue;
                }

                // for each cell perform the near neighbour search
                for (int i = 0; i < cellList[x][y][z].size(); i++) {
                    int currentIndex = cellList[x][y][z][i];
                    
                    // calculate the number of neighbours in all surrounding cells
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                int newX = x + dx;
                                int newY = y + dy;
                                int newZ = z + dz;

                                // calculate number of neighbours for all surrounding cells
                                if (newX >= 0 && newX < cellList.size() && newY >= 0 && newY < cellList[x].size() && newZ >= 0 && newZ < cellList[x][y].size()) {
                                    calculateNeighboursPerCell(cellList[newX][newY][newZ], currentIndex, rc2);
                                }
                            }
                        }
                    }

                }

                counter++;
            }
        }
    }

    // Use MPI_Reduce to sum all elements in the vector into the root vector
    std::vector<int> globalNumNeighbours(numParticles, 0);
    MPI_Reduce(&numNeighbours[0], &globalNumNeighbours[0], numParticles, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // if the current process is the root, set the global vector to the root vector
    if (world_rank == 0) {
        numNeighbours = globalNumNeighbours;
    }
}


void System::calculateNumNeighboursOpenMPStatic(double rc) {
    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);

    // square rc so that it only has to be calculated once
    double rc2 = rc*rc;

    // loop over the particles
    for (int i=0;i<numParticles;i++) {
        // loop from i+1 to numParticles to reduce the number of calculations
        #pragma omp parallel for schedule(static)
        for (int j=i+1;j<numParticles;j++) {
            double r2 = 0.0;
            for (int k=0;k<3;k++) {
                r2 += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
            }
            if (r2<rc2) {
                #pragma omp atomic
                numNeighbours[i]++;
                #pragma omp atomic
                numNeighbours[j]++;
            }
        }
    }
}


void System::calculateNumNeighboursOpenMPDynamic(double rc) {
    // populate the vector and set it to zero for efficiency
    numNeighbours.resize(numParticles, 0);

    // square rc so that it only has to be calculated once
    double rc2 = rc*rc;

    // loop over the particles
    for (int i=0;i<numParticles;i++) {
        // loop from i+1 to numParticles to reduce the number of calculations
        #pragma omp parallel for schedule(dynamic)
        for (int j=i+1;j<numParticles;j++) {
            double r2 = 0.0;
            for (int k=0;k<3;k++) {
                r2 += (positions[i][k]-positions[j][k])*(positions[i][k]-positions[j][k]);
            }
            if (r2<rc2) {
                #pragma omp atomic
                numNeighbours[i]++;
                #pragma omp atomic
                numNeighbours[j]++;
            }
        }
    }
}


void System::calculateNeighboursPerCellOpenMP(std::vector <int> currentCell, int particleIndex, double rc2) {
    // loop over the particles in the current cell
    for (int j = 0; j < currentCell.size(); j++) {
        int neighbourIndex = currentCell[j];
        if (particleIndex != neighbourIndex) {
            double r2 = 0.0;
            for (int k = 0; k < 3; k++) {
                r2 += (positions[particleIndex][k] - positions[neighbourIndex][k]) * (positions[particleIndex][k] - positions[neighbourIndex][k]);
            }
            if (r2 < rc2) {
                // use an atomic operation to increment the number of neighbours
                #pragma omp atomic
                numNeighbours[particleIndex]++;
            }
        }
    }
}


void System::calculateNumNeighboursOpenMPCellList(double rc) {
    numNeighbours.resize(numParticles, 0);

    std::vector< std::vector< std::vector< std::vector<int> > > > cellList = createCellList(rc);

    double rc2 = rc * rc;

    // loop through the 3D cell list
    #pragma omp parallel for
    for (int x = 0; x < cellList.size(); x++) {
        for (int y = 0; y < cellList[x].size(); y++) {
            for (int z = 0; z < cellList[x][y].size(); z++) {
                
                // for each cell perform the near neighbour search
                for (int i = 0; i < cellList[x][y][z].size(); i++) {
                    int currentIndex = cellList[x][y][z][i];
                    
                    // calculate the number of neighbours in all surrounding cells
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                int newX = x + dx;
                                int newY = y + dy;
                                int newZ = z + dz;

                                // calculate number of neighbours for all surrounding cells
                                if (newX >= 0 && newX < cellList.size() && newY >= 0 && newY < cellList[x].size() && newZ >= 0 && newZ < cellList[x][y].size()) {
                                    calculateNeighboursPerCellOpenMP(cellList[newX][newY][newZ], currentIndex, rc2);
                                }
                            }
                        }
                    }

                }

            }
        }
    }
}

