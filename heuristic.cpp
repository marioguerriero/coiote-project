#include <iostream>
#include <random>
#include <thread>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <sys/time.h>
#include "heuristic.h"

#include <math.h>
#include <cstdlib>
#include <ctime>
#include <unordered_set>

using namespace std;

Data Heuristic::problem = Data();
int**** Heuristic::solution;
int Heuristic::demand;
int Heuristic::nCells;
int Heuristic::nCustomerTypes;
int Heuristic::nTimeSteps;
double Heuristic::execTime = 0;
double Heuristic::execTimeStart = 0;

Heuristic::Heuristic(string path){
    this->hasSolution = false;
    string line;
    string word;

    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        cin.get();
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word;
    this->nCells = atoi(word.c_str());
    iss >> word;
    this->nTimeSteps = atoi(word.c_str());
    iss >> word;
    this->nCustomerTypes = atoi(word.c_str());

    // Memory allocation
    solution = new int***[nCells];
    problem.costs = new double***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.costs[i] = new double**[nCells];
        solution[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            problem.costs[i][j] = new double*[nCustomerTypes];
            solution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                problem.costs[i][j][m] = new double[nTimeSteps];
                solution[i][j][m] = new int[nTimeSteps];
            }
        }
    }
    problem.n = new int[nCustomerTypes];
    problem.activities = new int[nCells];
    problem.usersCell = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.usersCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            problem.usersCell[i][m] = new int[nTimeSteps];
        }
    }

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issN(line);
    for (int m = 0; m < nCustomerTypes; m++) {
        issN >> word;
        problem.n[m] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);// linea con m e t
            for (int i = 0; i < nCells; i++) {
                getline(iffN, line);// linea della matrice c_{ij} per t ed m fissati
                istringstream issC(line);
                for (int j = 0; j < nCells; j++) {
                    issC >> word;
                    problem.costs[i][j][m][t] = atoi(word.c_str());
                }
            }
        }
    }

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issA(line);
    for (int i = 0; i < nCells; i++) {
        issA >> word;
        problem.activities[i] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);
            getline(iffN, line);
            std::replace(line.begin(), line.end(), ';', ' ');
            istringstream issU(line);
            for (int i = 0; i < nCells; i++) {
                issU >> word;
                problem.usersCell[i][m][t] = atoi(word.c_str());
            }
        }
    }

    iffN.close();
    iss.str("");
}

vector<vector<int *>> optimalCombinations((unsigned long) Heuristic::nCells, vector<int*>(10));
int counter = 0;

double Heuristic::solveFast(vector<double>& stat, int timeLimit, float pRep, float pMut, float pCross) {
    //initializing variables of the problem
    execTimeStart = clock();

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;


    // Obtain all "certain moves" for all destinations

    // At each index we have a vector containing the optimal moves for the corresponding destination
    for(int j = 0; j < nCells; j++) {
        int demand = problem.activities[j];
        // Find the optimal combinations for the j-th cell
        vector<int *> combinations = getUsersCombination(demand);
        while (combinations.size() == 0)
            combinations = getUsersCombination(++demand);

        // Add the found combinations to the optimalCombinations vector's corresponding entry
        optimalCombinations.push_back(combinations);
    }

    //// Simulated Annealing

    // Initialization
    int ****currSolution = new int***[nCells];
    int ****bestSolution = new int***[nCells];
    for (int i = 0; i < nCells; i++) {
        currSolution[i] = new int**[nCells];
        bestSolution[i] = new int**[nCells];
        for (int j = 0; j < nCells; j++) {
            currSolution[i][j] = new int*[nCustomerTypes];
            bestSolution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < nCustomerTypes; m++) {
                currSolution[i][j][m] = new int[nTimeSteps];
                bestSolution[i][j][m] = new int[nTimeSteps];
                for(int t = 0; t < nTimeSteps; t++) {
                    currSolution[i][j][m][t] = 0;
                    bestSolution[i][j][m][t] = 0;
                }
            }
        }
    }

    int ***currUsers = new int**[nCells];
    int ***bestUsers = new int**[nCells];
    for (int i = 0; i < nCells; i++) {
        currUsers[i] = new int*[nCustomerTypes];
        bestUsers[i] = new int*[nCustomerTypes];
        for (int m = 0; m < nCustomerTypes; m++) {
            currUsers[i][m] = new int[nTimeSteps];
            bestUsers[i][m] = new int[nTimeSteps];
            for(int t = 0; t < nTimeSteps; t++) {
                currUsers[i][m][t] = problem.usersCell[i][m][t];
                bestUsers[i][m][t] = problem.usersCell[i][m][t];
            }
        }
    }

    //////// Greedy
    // We use a greedy heuristic to initialize the SA's initial population
    greedy(currSolution, currUsers);

    copySolution(bestSolution, currSolution);
    copyUsers(bestUsers, currUsers);

    int iteration = 0;

    // Simulated annealing iterations
    while((clock() - execTimeStart) / CLOCKS_PER_SEC < 4.8) { // stopping condition is 5 seconds
        neighboor(currSolution, currUsers); // Neighboor generation

        iteration++;

        double currentObj = objective(currSolution);
        double bestObj = objective(bestSolution);

        if(currentObj < bestObj) {
            copySolution(bestSolution, currSolution);
            copyUsers(bestUsers, currUsers);
        }
        else {
            double expP = exp((float) (-(currentObj-bestObj)/temperature(iteration)));
            double rndP = ((double) rand() / (1)) + 1;
            if(rndP < expP) {
                copySolution(bestSolution, currSolution);
                copyUsers(bestUsers, currUsers);
            }
        }
    }

    double objFun = objective(bestSolution);
    copySolution(solution, bestSolution);
    copyUsers(problem.usersCell, bestUsers);

    execTime = (clock() - execTimeStart) / CLOCKS_PER_SEC;
    stat.push_back(execTime);
    stat.push_back(objFun);
    hasSolution = true;
    return objFun;
}

void Heuristic::neighboor(int ****sol, int ***users) {
    /**
     * Pick two random m values m1 and m2. If they are one the multiple of the other than
     * pick a random i,j,t and swap a random amount of users with those of another cell
     */

    // Pick two random m values
    int m1 = (rand() % nCustomerTypes);
    int m2 = (rand() % nCustomerTypes);
    while(m2 == m1) m2 = (rand() % nCustomerTypes);

    // If they are one the multiple of the other continue
    if(problem.n[m1] % problem.n[m2] == 0 || problem.n[m2] % problem.n[m1] == 0) {
        // Generate a random i which will be our source cell
        int sourceI = (rand() % nCells);
        int sourceJ = (rand() % nCells);
        int sourceT = (rand() % nTimeSteps);

        // Choose source and destination types
        /** Generally the sourceM is the one with the higher problem.n.
         * Making this choice is pointless. It's just to simplify the next
         * steps of the algorithm.
         */
        int sourceM, destM;
        if(problem.n[m1] > problem.n[m2]) {
            sourceM = m1;
            destM = m2;
        }
        else {
            sourceM = m2;
            destM = m1;
        }

        // Find out how many users of type sourceM we want to move
        if(sol[sourceI][sourceJ][sourceM][sourceT] == 0) return;
        int usersToMoveSrc = 1 + (rand() % (sol[sourceI][sourceJ][sourceM][sourceT]));

        // Find out how many users of type destM we should move to fit the usersToMove value
        int usersToMoveDest = (usersToMoveSrc * problem.n[sourceM]) / problem.n[destM];

        // Chose the first destination capable of containing the amount of users we want to move
        int destI = -1, destJ = -1, destT = -1;
        bool found = false;
        for(int i = 0; i < nCells && !found; i++) {
            for(int j = 0; j < nCells && !found; j++) {
                for(int t = 0; t < nTimeSteps && !found; t++) {
                    if(i == sourceI && j == sourceJ && t == sourceT) continue;
                    if(sol[i][j][destM][t] >= usersToMoveSrc) {
                        destI = i, destJ = j, destT = t;
                        found = true;
                    }
                }
            }
        }
        if(!found) return;

        // Swap!
        sol[sourceI][sourceJ][sourceM][sourceT] -= usersToMoveSrc;
        sol[destI][destJ][sourceM][destT] += usersToMoveSrc;
        sol[destM][destJ][destM][destT] -= usersToMoveSrc;
        sol[sourceI][sourceJ][destM][sourceT] += usersToMoveSrc;

        users[sourceI][sourceM][sourceT] -= usersToMoveSrc;
        users[destI][destM][destT] -= usersToMoveDest;
    }
}

/**
 * This method build a feasible optimum solution with a good objective function
 * value which will be used as a starting point for our simulated annealing algorithm.
 * @param sol
 */
void Heuristic::greedy(int ****sol, int ***users) {
    // Initialize data structures
    //activities
    int *act = new int[nCells];
    for (int i = 0; i < nCells; i++) act[i] = problem.activities[i];

    int demand = 0;
    for(int i = 0; i < nCells; i++) demand += problem.activities[i];

    // Perform the greedy heuristic

    while(demand > 0) {
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;

        for (int i = 0; i < nCells; i++) {
            for (int j = 0; j < nCells; j++) {
                if(i == j || act[j] == 0) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    if(act[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(users[i][m][t]==0) continue;

                        if(minCost > ((double)problem.costs[i][j][m][t]/(double)problem.n[m])
                           || (minCost == ((double)problem.costs[i][j][m][t]/(double)problem.n[m]) && m > mt)) {
                            minCost = ((double)problem.costs[i][j][m][t]/(double)problem.n[m]);
                            it = i; mt = m; tt = t; jt = j;
                        }
                    }
                }
            }
        }

        // update solution1s
        int maxTasks = problem.n[mt] * users[it][mt][tt];

        if(maxTasks > act[jt]) {
            int movedUsers = act[jt] / problem.n[mt];
            sol[it][jt][mt][tt] += movedUsers;
            act[jt] -= movedUsers * problem.n[mt];
            users[it][mt][tt] -= movedUsers;
            demand -= movedUsers * problem.n[mt];
        }
        else {
            sol[it][jt][mt][tt] += users[it][mt][tt];
            act[jt] -= maxTasks;
            users[it][mt][tt] = 0;
            demand -= maxTasks;
        }
    }

    swapsolutions(sol, users);
}

/**
 * goal of this method is the improvement of the found solution by swapping solutions (user moves)
 *
 * it is divided in two parts:
 * 1) search of the best swap (the one which improves the more obj func)
 * 2) apply swap moving as much users as we can
 *
 * a particular solutions is represented by a move from i to j of a certain number of m-t users.
 * in order to find which is the best swap to do we calc i->j + i2->j2 cost (actual solution) and
 * combined cost i->j2 + i2->j, if the difference between those costs is convenient we do the swap
 *
 * one of the main strategies which improves the effectivness of the algorithm is the possibility to
 * swap solutions with users that still are in their origin cells
 *
 * NB: we can swap two same-type users or also two different-type users if there is a relation between them
 * (for example I can swap one 6 tasks user with two 3 tasks users). this kind of control is done in part (1)
 *
 */
void Heuristic::swapsolutions(int**** solutions, int*** users) {
    /*
     * part 1
     */
    //we want to maximize the cost difference (maxdiff) between the actual solution (>) and the future one (<)
    double maxdiff = 0;
    //loop until max costs difference is > 0
    do {
        maxdiff = 0;
        //temp variables needed below to apply founded solutions
        int idiff=-1, jdiff=-1, mdiff=-1, tdiff=-1, i2diff=-1, j2diff=-1, t2diff=-1, m2diff=-1;
        //flag used to keep track of which "type" of swap we have to do (same type users or different ones)
        bool swaptype; // 0 = same type.             1 = different type
        //flag used to differenciate "solutions to solutions" and "solutions to users" swap
        bool userswap; // 0 = s <-> s.              1 = s <-> u

        //start search of the max cost difference
        for (int j = 0; j < nCells; j++) {
            if (problem.activities[j] == 0) continue;
            for (int i = 0; i < nCells; i++) {
                if (i == j) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    for (int t = 0; t < nTimeSteps; t++) {
                        if (solutions[i][j][m][t] == 0) continue;

                        // i -> j cost
                        int ijcost = (int)problem.costs[i][j][m][t];

                        for(int m2 = 0; m2 < nCustomerTypes; m2++) {
                            //those are to check if it is possible to swap two different type users
                            int quotient = problem.n[m2] / problem.n[m];
                            int residual = problem.n[m2] % problem.n[m];
                            //if users are of the same type (or a combinable one) residual is 0
                            if (residual != 0) continue;

                            if (quotient == 1) {
                                //quotient == 1 -> m2 == m: swap one to one.
                                //restart search of the max difference
                                for (int j2 = 0; j2 < nCells; j2++) {
                                    if (problem.activities[j2] == 0) continue;
                                    // i -> j2 cost
                                    int ij2cost = (int)problem.costs[i][j2][m][t];
                                    for (int i2 = 0; i2 < nCells; i2++) {
                                        for (int t2 = 0; t2 < nTimeSteps; t2++) {
                                            //let's differenciate search on solutions and users
                                            if (solutions[i2][j2][m][t2] != 0) {
                                                // i2 -> j2 cost
                                                int i2j2cost = (int) problem.costs[i2][j2][m][t2];
                                                // i2 -> j cost
                                                int i2jcost = (int) problem.costs[i2][j][m][t2];

                                                double tempdiff = (ijcost + i2j2cost) - (i2jcost + ij2cost);
                                                if (maxdiff < tempdiff) {
                                                    idiff = i; jdiff = j; mdiff = m; tdiff = t;
                                                    i2diff = i2; t2diff = t2; j2diff = j2;
                                                    maxdiff = tempdiff;
                                                    swaptype = false; //flag false = same type swap
                                                    userswap = false;
                                                }
                                            }
                                            if(users[i2][m][t2] != 0) {
                                                int i2jcost = (int) problem.costs[i2][j][m][t2];
                                                double tempdiff = ijcost - i2jcost;
                                                if (maxdiff < tempdiff) {
                                                    idiff = i; jdiff = j; mdiff = m; tdiff = t;
                                                    i2diff = i2; t2diff = t2;
                                                    maxdiff = tempdiff;
                                                    swaptype = false; //flag false = same type swap
                                                    userswap = true;
                                                }
                                            }
                                        }
                                    }
                                }
                            } // else I need to check if there are enough users into i cell
                            else if(solutions[i][j][m][t] >= quotient) {
                                for (int j2 = 0; j2 < nCells; j2++) {
                                    if (problem.activities[j2] == 0) continue;
                                    int ij2cost = (int)problem.costs[i][j2][m][t];
                                    for (int i2 = 0; i2 < nCells; i2++) {
                                        for (int t2 = 0; t2 < nTimeSteps; t2++) {
                                            if (solutions[i2][j2][m2][t2] != 0) {
                                                int i2j2cost = (int) problem.costs[i2][j2][m2][t2];
                                                int i2jcost = (int) problem.costs[i2][j][m2][t2];

                                                double tempdiff = (quotient * ijcost + i2j2cost) - (i2jcost + quotient * ij2cost);
                                                if (maxdiff < tempdiff) {
                                                    idiff = i; jdiff = j; mdiff = m; tdiff = t;
                                                    i2diff = i2; t2diff = t2; j2diff = j2; m2diff = m2;
                                                    maxdiff = tempdiff;
                                                    swaptype = true; //flag true = different type swap
                                                    userswap = false;
                                                }
                                            }
                                            if(users[i2][m2][t2] != 0) {
                                                int i2jcost = (int) problem.costs[i2][j][m2][t2];
                                                double tempdiff = quotient * ijcost - i2jcost;
                                                if (maxdiff < tempdiff) {
                                                    idiff = i; jdiff = j; mdiff = m; tdiff = t;
                                                    i2diff = i2; t2diff = t2; m2diff = m2;
                                                    maxdiff = tempdiff;
                                                    swaptype = true; //flag true = different type swap
                                                    userswap = true;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //break if the swap will not improve solution anymore
        if(maxdiff <= 0) break;

        /*
         * part 2
         */
        //while until there are enough users to move
        if(swaptype) {
            if(userswap) {
                while (solutions[idiff][jdiff][mdiff][tdiff] >= (problem.n[m2diff] / problem.n[mdiff]) &&
                       users[i2diff][m2diff][t2diff] != 0) {
                    solutions[idiff][jdiff][mdiff][tdiff] -= problem.n[m2diff] / problem.n[mdiff];
                    users[idiff][mdiff][tdiff] += problem.n[m2diff] / problem.n[mdiff];
                    solutions[i2diff][jdiff][m2diff][t2diff]++;
                    users[i2diff][m2diff][t2diff]--;
                }
            } else {
                while (solutions[idiff][jdiff][mdiff][tdiff] >= (problem.n[m2diff] / problem.n[mdiff]) &&
                       solutions[i2diff][j2diff][m2diff][t2diff] != 0) {
                    solutions[idiff][jdiff][mdiff][tdiff] -= problem.n[m2diff] / problem.n[mdiff];
                    solutions[i2diff][j2diff][m2diff][t2diff]--;
                    solutions[idiff][j2diff][mdiff][tdiff] += problem.n[m2diff] / problem.n[mdiff];
                    solutions[i2diff][jdiff][m2diff][t2diff]++;
                }
            }
        } else {
            if(userswap) {
                while (solutions[idiff][jdiff][mdiff][tdiff] != 0 && users[i2diff][mdiff][t2diff] != 0) {
                    solutions[idiff][jdiff][mdiff][tdiff]--;
                    users[i2diff][mdiff][t2diff]--;
                    solutions[i2diff][jdiff][mdiff][t2diff]++;
                    users[idiff][mdiff][tdiff]++;
                }
            } else {
                while (solutions[idiff][jdiff][mdiff][tdiff] != 0 && solutions[i2diff][j2diff][mdiff][t2diff] != 0) {
                    solutions[idiff][jdiff][mdiff][tdiff]--;
                    solutions[i2diff][j2diff][mdiff][t2diff]--;
                    solutions[idiff][j2diff][mdiff][tdiff]++;
                    solutions[i2diff][jdiff][mdiff][t2diff]++;
                }
            }
        }

    } while(maxdiff > 0);

}

double Heuristic::temperature(double input) {
    return 0.999*input;
}

double Heuristic::objective(int ****sol) {
    double obj = 0.0;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if(sol[i][j][m][t] > 0) obj += sol[i][j][m][t] * problem.costs[i][j][m][t];
    return obj;
}

void Heuristic::writeKPI(string path, string instanceName, vector<double> stat){
    if (!hasSolution)
        return;

    ofstream fileO(path, ios::app);
    if(!fileO.is_open())
        return;

    fileO << instanceName << ";" << stat[0] << ";" << stat[1];
    for(int i=2; i<stat.size(); i++)
        fileO <<  ";" << stat[i];
    fileO << endl;

    fileO.close();

}

void Heuristic::writeSolution(string path) {
    if (!hasSolution)
        return;

    ofstream fileO(path);
    if(!fileO.is_open())
        return;

    fileO << this->nCells << "; " << this->nTimeSteps << "; " << this->nCustomerTypes << endl;
    for (int m = 0; m < this->nCustomerTypes; m++)
        for (int t = 0; t < this->nTimeSteps; t++)
            for (int i = 0; i < this->nCells; i++)
                for (int j = 0; j < this->nCells; j++)
                    if (solution[i][j][m][t] > 0)
                        fileO << i << ";" << j << ";" << m << ";" << t << ";" << solution[i][j][m][t] << endl;

    fileO.close();
}

eFeasibleState Heuristic::isFeasible(string path) {

    string line;
    string word;
    int nCellsN;
    int nTimeStepsN;
    int nCustomerTypesN;
    int i, j, m, t;


    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word; // nCells
    nCellsN = atoi(word.c_str());
    iss >> word; // nTimeSteps
    nTimeStepsN = atoi(word.c_str());
    iss >> word; // nCustomerTypes
    nCustomerTypesN = atoi(word.c_str());

    int**** solutionN = new int***[nCells];
    for (i = 0; i < nCellsN; i++) {
        solutionN[i] = new int**[nCells];
        for (j = 0; j < nCellsN; j++) {
            solutionN[i][j] = new int*[nCustomerTypes];
            for (m = 0; m < nCustomerTypesN; m++) {
                solutionN[i][j][m] = new int[nTimeSteps];
                for ( t = 0; t < nTimeStepsN; t++) {
                    solutionN[i][j][m][t] = 0;
                }
            }
        }
    }

    while (getline(iffN, line)) {
        std::replace(line.begin(), line.end(), ';', ' ');
        istringstream iss(line);
        iss >> word; // i
        i = atoi(word.c_str());
        iss >> word; // j
        j = atoi(word.c_str());
        iss >> word; // m
        m = atoi(word.c_str());
        iss >> word; // t
        t = atoi(word.c_str());
        iss >> word; // value
        solutionN[i][j][m][t] = atoi(word.c_str());
    }

    // Demand
    bool feasible = true;
    int expr;
    for (int i = 0; i < nCells; i++) {
        expr = 0;
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    expr += problem.n[m] * solutionN[j][i][m][t];
        if (expr < problem.activities[i])
            feasible = false;
    }

    if (!feasible)
        return NOT_FEASIBLE_DEMAND;

    // Max Number of users
    for (int i = 0; i < nCells; i++)
        for (int m = 0; m < nCustomerTypes; m++)
            for (int t = 0; t < nTimeSteps; t++) {
                expr = 0;
                for (int j = 0; j < nCells; j++)
                    expr += solutionN[i][j][m][t];
                if (expr > problem.usersCell[i][m][t])
                    feasible = false;
            }

    if(!feasible)
        return NOT_FEASIBLE_USERS;

    return FEASIBLE;
}

void Heuristic::getStatSolution(vector<double>& stat) {
    if (!hasSolution)
        return;

    int* tipi = new int[nCustomerTypes];
    for (int m = 0; m < nCustomerTypes; m++)
        tipi[m] = 0;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int t = 0; t < nTimeSteps; t++)
                for (int m = 0; m < nCustomerTypes; m++)
                    if (solution[i][j][m][t] > 0)
                        tipi[m] += solution[i][j][m][t];
    for (int m = 0; m < nCustomerTypes; m++)
        stat.push_back(tipi[m]);

}

//TODO: order problem.n[] and combinations
vector<int*> Heuristic::getUsersCombination(int demand) {
    vector<int*> combinations;
    int* indexes = new int[nCustomerTypes];
    bool* initializeIndexes = new bool[nCustomerTypes];
    int* lastDemands = new int[nCustomerTypes];
    int* usersChosen = new int[nCustomerTypes];
    for(int m=0; m<nCustomerTypes; m++){
        indexes[m] = 0;
        initializeIndexes[m] = true;
        lastDemands[m] = 0;
        usersChosen[m] = 0;
    }
    int lastUserType = nCustomerTypes-1;
    while(true){
        if(indexes[lastUserType] < 0){
            if(lastUserType == nCustomerTypes-1) break;
            indexes[lastUserType] = 0;
            lastUserType++;
            continue;
        }
        if(initializeIndexes[lastUserType]){
            if(lastUserType == nCustomerTypes-1) indexes[lastUserType] = demand/problem.n[lastUserType];
            else indexes[lastUserType] = lastDemands[lastUserType+1]/problem.n[lastUserType];
            initializeIndexes[lastUserType] = false;
        }
        usersChosen[lastUserType] = indexes[lastUserType];
        if(lastUserType == nCustomerTypes-1) lastDemands[lastUserType] = demand - (usersChosen[lastUserType] * problem.n[lastUserType]);
        else lastDemands[lastUserType] = lastDemands[lastUserType+1] - (usersChosen[lastUserType] * problem.n[lastUserType]);
        if(lastUserType == 0){
            if(lastDemands[lastUserType] == 0) {
                int *combination = new int[nCustomerTypes];
                for(int i = 0; i < nCustomerTypes; i++)
                    combination[i] = usersChosen[i];
                combinations.push_back(combination);
                //if (demand == 22) cout << usersChosen[2] << " " << usersChosen[1] << " " << usersChosen[0] << endl;
            }
            lastUserType++;
            continue;
        }
        indexes[lastUserType]--;
        lastUserType--;
        initializeIndexes[lastUserType] = true;
    }

    return combinations;
}

void Heuristic::copySolution(int ****dest, int ****src) {
    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    dest[i][j][m][t] = src[i][j][m][t];
}

void Heuristic::copyUsers(int ***dest, int ***src) {
    for (int i = 0; i < nCells; i++)
        for (int m = 0; m < nCustomerTypes; m++)
            for (int t = 0; t < nTimeSteps; t++)
                dest[i][m][t] = src[i][m][t];
}

/*
vector<int *> Heuristic::greedy(vector<vector<int *>> exactCombinations) {
    // Make local copies of the problem variables
    int demand = 0;
    for(int i = 0; i < nCells; i++) {
        demand += problem.activities[i];
    }

    int ****solution = new int***[nCells];
    for (int i = 0; i < nCells; i++) {
        solution[i] = new int**[nCells];
        for (int j = 0; j < nCells; j++) {
            solution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < nCustomerTypes; m++) {
                solution[i][j][m] = new int[nTimeSteps];
                for (int t = 0; t < nTimeSteps; t++) {
                    solution[i][j][m][t] = 0;
                }
            }
        }
    }

    int ***usersCell = new int**[nCells];
    for (int i = 0; i < nCells; i++) {
        usersCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < nCustomerTypes; m++) {
            usersCell[i][m] = new int[nTimeSteps];
            for (int t = 0; t < nTimeSteps; t++) {
                usersCell[i][m][t] = problem.usersCell[i][m][t];
            }
        }
    }

    int *activities = new int[nCells];
    for (int i = 0; i < nCells; i++) {
        activities[i] = problem.activities[i];
    }

    // Perform only those moves that must be done in a certain way
    std::unordered_set<int> tabooList;

    for(int x = 0; x < nCells; x++) {
        // Find the j with the minimum number of activities
        int minActivities = activities[0];
        int minActivitiesJ = -1;

        for (int j = 0; j < nCells; j++) {
            if(!activities[j]) continue;
            if(tabooList.count(j)) continue;
            int currAct = activities[j];
            if (minActivities > currAct) {
                minActivities = currAct;
                minActivitiesJ = j;
            }
        }

        tabooList.insert(minActivitiesJ);
        if(minActivitiesJ == -1) continue;

        // Perform the moves for the found j
        // Find the minimum combination
        vector<int*> moves = exactCombinations[minActivitiesJ];
        int *minUsersToMove = new int[nCustomerTypes];
        for(int m = 0; m < nCustomerTypes; m++) {
            minUsersToMove[m] = moves[0][m];
        }

        for(int i = 0; i < moves.size(); i++) {
            int *currMove = moves[i];

            for(int m = 0; m < nCustomerTypes; m++) {
                int currMovedUsers = currMove[m];
                if(currMovedUsers < minUsersToMove[m])
                    minUsersToMove[m] = currMovedUsers;
            }
        }

        // Perform the minimum move
        for(int m = 0; m < nCustomerTypes; m++) {
            int users = minUsersToMove[m];
            if(users == 0) continue;

            // Look for a proper destination
            while(users > 0) {
                double minCost = 1e10;
                int minI = -1, minT = -1;
                for(int i = 0; i < nCells; i++) {
                    if(i == minActivitiesJ) continue;
                    for(int t = 0; t < nTimeSteps; t++) {
                        if(!usersCell[i][m][t]) continue;

                        double currCost = problem.costs[i][minActivitiesJ][m][t];
                        if(minCost > currCost) {
                            minCost = currCost;
                            minI = i; minT = t;
                        }
                    }
                }

                // Make the move
                if(users < usersCell[minI][m][minT]) {
                    usersCell[minI][m][minT] -= users;
                    activities[minActivitiesJ] -= users * problem.n[m];
                    demand -= users * problem.n[m];
                    solution[minI][minActivitiesJ][m][minT] += users;
                    users = 0;
                }
                else {
                    users -= usersCell[minI][m][minT];
                    activities[minActivitiesJ] -= usersCell[minI][m][minT] * problem.n[m];
                    demand -= usersCell[minI][m][minT] * problem.n[m];
                    solution[minI][minActivitiesJ][m][minT] += usersCell[minI][m][minT];
                    usersCell[minI][m][minT] = 0;
                }

                if(activities[minActivitiesJ] < 0) {
                    demand += -activities[minActivitiesJ];
                    activities[minActivitiesJ] = 0;
                }
            }
        }
    }

    // Perform the remaining moves

    // Calculate the available users
    int availableUsers[nCustomerTypes];
    for(int m = 0; m < nCustomerTypes; m++) availableUsers[m] = 0;
    for(int i = 0; i < nCells; i++)
        for (int m = 0; m < nCustomerTypes; m++)
            for (int t = 0; t < nTimeSteps; t++)
                availableUsers[m] += usersCell[i][m][t];


    // Recalculate combinations
    vector<vector<int*>> optComb;
    for(int j = 0; j < nCells; j++) {
        int dem = activities[j];
        // Find the optimal combinations for the j-th cell
        vector<int *> combinations = getUsersCombination(dem);
        while (combinations.size() == 0)
            combinations = getUsersCombination(++dem);

        // Add the found combinations to the optimalCombinations vector's corresponding entry
        optComb.push_back(combinations);
    }

    // Actually perform the moves
    while(demand > 0) {
        // Find the j with the minimum number of activities
        int minActivities = activities[0];
        int minActivitiesJ = 0;

        for (int j = 0; j < nCells; j++) {
            if(!activities[j]) continue;
            int currAct = activities[j];
            if (minActivities > currAct) {
                minActivities = currAct;
                minActivitiesJ = j;
            }
        }

        // Perform the moves for the found j
        // Find the best combination
        vector<int*> moves = optComb[minActivitiesJ];
        int *move = new int[nCustomerTypes];
        memset(move, (int) 1e5, nCustomerTypes*sizeof(int));
        int movingUsers = 0;
        for(int m = 0; m < nCustomerTypes; m++) movingUsers += move[m];

        // Find the move which moves the minimum number of users
        for(int i = 0; i < moves.size(); i++) {
            int *currMove = moves[i];

            bool enoughUsers = true;
            for(int m = 0; m < nCustomerTypes; m++) {
                enoughUsers = (availableUsers[m] >= currMove[m]);
                if(!enoughUsers) break;
            }
            if(!enoughUsers) continue;

            int currMovingUsers = 0;
            for(int m = 0; m < nCustomerTypes; m++) {
                currMovingUsers += currMove[m];
            }

            if(currMovingUsers < movingUsers) {
                movingUsers = currMovingUsers;
                move = currMove;
            }
        }

        // Perform the minimum move
        for(int m = 0; m < nCustomerTypes; m++) {
            int users = move[m];
            if(users == 0) continue;

            // Look for a proper destination
            while(users != 0) {
                if(users <= 0)
                    break;

                double minCost = 1e10;
                int minI = -1, minT = -1;
                for(int i = 0; i < nCells; i++) {
                    if(i == minActivitiesJ) continue;
                    for(int t = 0; t < nTimeSteps; t++) {
                        if(!usersCell[i][m][t]) continue;

                        double currCost = problem.costs[i][minActivitiesJ][m][t];
                        if(minCost > currCost) {
                            minCost = currCost;
                            minI = i; minT = t;
                        }
                    }
                }

                // Make the move
                if(users < usersCell[minI][m][minT]) {
                    usersCell[minI][m][minT] -= users;
                    activities[minActivitiesJ] -= users * problem.n[m];
                    demand -= users * problem.n[m];
                    solution[minI][minActivitiesJ][m][minT] += users;
                    availableUsers[m] -= users;
                    users = 0;
                }
                else {
                    users -= usersCell[minI][m][minT];
                    activities[minActivitiesJ] -= usersCell[minI][m][minT] * problem.n[m];
                    demand -= usersCell[minI][m][minT] * problem.n[m];
                    solution[minI][minActivitiesJ][m][minT] += usersCell[minI][m][minT];
                    availableUsers[m] -= usersCell[minI][m][minT];
                    usersCell[minI][m][minT] = 0;
                }

                if(activities[minActivitiesJ] < 0) {
                    demand += -activities[minActivitiesJ];
                    activities[minActivitiesJ] = 0;
                }
            }
        }
    }

    // Free memory
    delete solution;
    delete usersCell;
    delete activities;

    double obj = 0.0;
    for(int i = 0; i < nCells; i++)
        for(int j = 0; j < nCells; j++)
            for(int m = 0; m < nCustomerTypes; m++)
                for(int t = 0; t < nTimeSteps; t++)
                    obj += solution[i][j][m][t] * problem.costs[i][j][m][t];

    return vector<int *>();
}*/