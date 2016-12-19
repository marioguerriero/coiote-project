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
    double objFun = 0;
    execTimeStart = clock();
    timeval time;
    gettimeofday(&time, NULL);
    execTimeStart = time.tv_sec;

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
    for (int i = 0; i < this->nCells; i++) {
        currSolution[i] = new int**[nCells];
        bestSolution[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            currSolution[i][j] = new int*[nCustomerTypes];
            bestSolution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                currSolution[i][j][m] = new int[nTimeSteps];
                bestSolution[i][j][m] = new int[nTimeSteps];
            }
        }
    }

    //////// Greedy
    // We use a greedy heuristic to initialize the SA's initial population
    greedy(currSolution);

    copySolution(bestSolution, currSolution);

    int iteration = 0;

    // Simulated annealing iterations
    while((clock() - execTimeStart) / CLOCKS_PER_SEC < 5) { // stopping condition is 5 seconds
        neighboor(currSolution); // Neighboor generation

        iteration++;

        double currentObj = objective(currSolution);
        double bestObj = objective(bestSolution);

        if(currentObj < bestObj) {
            copySolution(bestSolution, currSolution);
        }
        else {
            double expP = exp((float) (-(currentObj-bestObj)/temperature(iteration)));
            double rndP = ((double) rand() / (1)) + 1;
            if(rndP < expP) {
                copySolution(bestSolution, currSolution);
            }
        }
    }

    execTime = (clock() - execTimeStart) / CLOCKS_PER_SEC;
    stat.push_back(execTime);
    stat.push_back(objFun);
    hasSolution = true;
    return objFun;
}

void Heuristic::neighboor(int ****sol) {
    /**
     * Generate a random x and a random (i,j,m,t) and swap x users from
     * (i,j,m,t) to the first fitting cell eg if sol[i][j][m][t] is equal to 4
     * then find a j where problem.activities[j] is at least equal to
     * problem.n[m] * sol[i][j][m][t]. Now move x users (of any type) from
     * the old (i,j,m,t) to the new one.
     */

    //// Initialization

    // Generate (i,j,m,t)
    int i = (rand() % (int)(nCells + 1));
    int j = (rand() % (int)(nCells + 1));
    int m = (rand() % (int)(nCustomerTypes + 1));
    int t = (rand() % (int)(nTimeSteps + 1));

    // Generate x
    int x = (rand() % (int)(problem.usersCell[i][m][t] + 1));

    //// Swap

    // Count the users involved in the swap
    int usersToMove = sol[i][j][m][t];

    // Find the destination cell
    int destination = -1;
    for(int c = 0; c < nCells; c++) {
        if(c == j) continue;
        if(problem.activities[c] >= usersToMove * problem.n[m]) {
            destination = c;
            break;
        }
    }
    // If not destination could be found return
    if(destination == -1)
        return;

    // Find the best combinations both for j and for destination
    vector<int*> combJ = getUsersCombination(problem.activities[j]);
    vector<int*> combDest = getUsersCombination(problem.activities[destination]);


    // If there are two combinations where c1[m] is equal c2[m], then swap
    for(int k = 0; k < combJ.size(); k++) {
        int *currJComb = combJ[k];
        for(int w = 0; w < combDest.size(); w++) {
            int *currDestComb = combDest[w];
            if(currJComb[m] == x && currJComb[m] == currDestComb[m]) {

            }
        }
    }

}

void Heuristic::greedy(int ****sol) {
    /**
     * Find a set of optimum combinations and move users according to it in order
     * to find a feasible solution. That feasible solution will be our starting point
     * for the whole algorithm.
     */

    // Find a set of optimum combinations
    vector<vector<int*>> optCombs;
    for(int j = 0; j < nCells; j++) {
        int demand = problem.activities[j];
        optCombs.push_back(getUsersCombination(demand));
    }
}

/*double Heuristic::objective(vector<int *> population) {
    double obj = 0;

    // Initialize data structures used to analyze the objective value

    // Check how many available users are there
    int *availableUsers = new int[nCustomerTypes];
    memset(availableUsers, 0, nCustomerTypes*sizeof(int));
    for(int i = 0; i < nCells; i++)
        for(int t = 0; t < nTimeSteps; t++)
            for(int m = 0; m < nCustomerTypes; m++)
                availableUsers[m] += problem.usersCell[i][m][t];

    // solutions
    int ****sol = new int***[nCells];
    for (int i = 0; i < nCells; i++) {
        sol[i] = new int **[nCells];
        for (int j = 0; j < nCells; j++) {
            sol[i][j] = new int *[nCustomerTypes];
            for (int m = 0; m < nCustomerTypes; m++) {
                sol[i][j][m] = new int[nTimeSteps];
                for (int t = 0; t < nTimeSteps; t++)
                    sol[i][j][m][t] = 0;
            }
        }
    }

    // usersCell
    int ***uCell = new int**[nCells];
    for (int i = 0; i < nCells; i++) {
        uCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < nCustomerTypes; m++) {
            uCell[i][m] = new int[nTimeSteps];
            for (int t = 0; t < nTimeSteps; t++)
                uCell[i][m][t] = problem.usersCell[i][m][t];
        }
    }

    // activities
    int *act = new int[nCells];
    for (int i = 0; i < nCells; i++)
        act[i] = problem.activities[i];

    // Move users for each found combination
    for(int j = 0; j < population.size(); j++) {
        int *currComb = population[j];

        for (int m = 0; m < nCustomerTypes; m++) {
            if (currComb[m] == 0) continue;

            double minCost;
            while (currComb[m] > 0) {
                if(availableUsers[m] <= 0)
                    return 0.0;

                minCost = 1e10;
                int it = -1;
                int tt = -1;
                for (int i = 0; i < nCells; i++) {
                    if (i == j) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if (uCell[i][m][t] <= 0) continue;
                        if(availableUsers[m] <= 0) continue;
                        if (minCost > (problem.costs[i][j][m][t])) {
                            minCost = (problem.costs[i][j][m][t]);
                            it = i;
                            tt = t;
                        }
                    }
                }

                if (uCell[it][m][tt] >= currComb[m]) {
                    demand -= currComb[m] * problem.n[m];
                    act[j] -= currComb[m] * problem.n[m];
                    sol[it][j][m][tt] += currComb[m];
                    uCell[it][m][tt] -= currComb[m];
                    availableUsers[m] -= currComb[m];
                    currComb[m] = 0;
                } else {
                    sol[it][j][m][tt] += uCell[it][m][tt];
                    act[j] -= uCell[it][m][tt] * problem.n[m];
                    demand -= uCell[it][m][tt] * problem.n[m];
                    currComb[m] -= uCell[it][m][tt];
                    availableUsers[m] -= uCell[it][m][tt];
                    uCell[it][m][tt] = 0;
                }

                if (act[j] < 0) {
                    demand += -act[j];
                    availableUsers[m] += -act[j];
                    act[j] = 0;
                }
            }
        }
    }

    // Calculate the objective function value
    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if(sol[i][j][m][t] > 0) obj += sol[i][j][m][t] * problem.costs[i][j][m][t];

    cout << "Objective function: " << obj << endl;

    return obj;
}*/

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