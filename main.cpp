#include "system.h"



void Log(const char* filename, unsigned int N, double L, int choice, int numProc, double timeTaken) {
    /**
    * Logs the given information to a file.
    * 
    * filename The name of the file to log the information to.
    * N The number of particles.
    * L The length of the cube containing the particles.
    * choice The choice value.
    * numProc The number of processors.
    * timeTaken The time taken.
    */
    
    std::ifstream inFile(filename, std::ios::in);
    std::string firstline;

    // Check if the file was opened successfully
    if (inFile.is_open()) {
        std::getline(inFile, firstline);
    }
    inFile.close();

    std::string header = "N,L,choice,numProc,numThread,timeTaken";

    if (firstline != header) {
        std::ofstream outFile(filename, std::ios::out);
        outFile << header << std::endl
        << N << ',' << L << ',' << choice << ',' << numProc << ',' << omp_get_max_threads() << ',' << timeTaken << std::endl;
        outFile.close();
    } else {
        std::ofstream outFile(filename, std::ios::app);
        outFile << N << ',' << L << ',' << choice << ',' << numProc << ',' << omp_get_max_threads() << ',' << timeTaken << std::endl;
        outFile.close();
    }
}


// entry point, takes in the filename and the choice of method and the cutoff radius

int main(int argc, char** argv) {
    /**
    * Entry point of the program.
    *
    * This program calculates the number of neighbors for each particle in a system.
    * It supports different calculation methods based on the choice provided as a command line argument.
    *
    * Usage: ./main <filename> <choice> <rc>
    *   - filename: Path to the input file containing system information.
    *   - choice: Calculation method choice (0-6).
    *   - rc: Cutoff radius for neighbor calculation.
    *
    * Compute Methods:
    *   0: Serial calculation
    *   1: Optimised serial calculation
    *   2: MPI calculation
    *   3: Optimised MPI calculation
    *   4: Balanced optimised MPI calculation
    *   5: OpenMP static calculation
    *   6: OpenMP dynamic calculation
    *
    *  argc The number of command line arguments.
    *  argv An array of command line arguments.
    *  returns 0 on success, non-zero on failure.
    */
    
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    char* filename = argv[1];
    int choice = atoi(argv[2]);
    double rc = atof(argv[3]);

    unsigned int N;
    double L;

    System system(filename, N, L);

    // measure the walltime
    double startTime = MPI_Wtime();

    // Call the function to calculate the number of neighbours based on the choice
    if (choice == 0) {
        system.calculateNumNeighboursSerial(rc);
    } else if (choice == 1) {
        system.calculateNumNeighboursSerialOptimised(rc);
    } else if (choice == 2) {
        system.calculateNumNeighboursMPI(rc, world_rank, world_size);
    } else if (choice == 3) {
        system.calculateNumNeighboursMPIOptimised(rc, world_rank, world_size);
    } else if (choice == 4) {
        system.calculateNumNeighboursMPIBalancedOptimised(rc, world_rank, world_size);
    } else if (choice == 5) {
        system.calculateNumNeighboursOpenMPStatic(rc);
    } else if (choice == 6) {
        system.calculateNumNeighboursOpenMPDynamic(rc);
    } else if (choice == 7) {
        system.calculateNumNeighboursSerialCellList(rc);
    } else if (choice == 8) {
        system.calculateNumNeighboursMPICellList(rc, world_rank, world_size);
    } else if (choice == 9) {
        system.calculateNumNeighboursOpenMPCellList(rc);
    }

    double endTime = MPI_Wtime();

    // root node logs the information
    if (world_rank == 0) {
        Log("log.csv", N, L, choice, world_size, endTime - startTime);

        for (int i = 0; i < N; i++) {
            std::cout << i << "  " << system.GetNumNeighbours()[i] << std::endl;
        }
    }

    MPI_Finalize();
}