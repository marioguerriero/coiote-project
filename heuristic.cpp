#include <iostream>
#include <random>
#include <thread>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "heuristic.h"

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


double Heuristic::solveFast(vector<double>& stat, int timeLimit, float pRep, float pMut, float pCross) {
    double objFun = 0;
    execTimeStart = clock();

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;
/*
    //calc total number of tasks
    int numTasks = 0;
    for(int i = 0; i < nCells; i++) {
        numTasks += problem.activities[i];
    }

    double singleTaskTime = timeLimit/numTasks;

    //j destination cell
    // upper level parallel heuristc
    for (int j = 0; j < nCells; j++) {
        // number of tasks to be done in the destination cell
        int demand = problem.activities[j];
        if(demand == 0) continue;

        PartialSolution **partialSolutions = new PartialSolution*[nCustomerTypes];

        for(int m = 0; m < nCustomerTypes; m++) {
            // time limit for each thread
            double timelimit = singleTaskTime * demand;

            //init second array (on m) of **partialSolutions
            partialSolutions[m] = new PartialSolution[nCells*nTimeSteps];
            // Begin thread of execution
            upperHeuristic(j, m, numTasks, timelimit, partialSolutions[m]);
        }

        PartialSolution *beamSolution = new PartialSolution[nTimeSteps*nCells*nCustomerTypes];
        int x = 0;
        for(int i = 0; i < nCustomerTypes; i++) {
            for (int k = 0; k < nCells; k++, x++) {
                beamSolution[x] = partialSolutions[i][k];
            }
        }

        int reminder = 0;
        while(demand > 0) {
            double minCost = 1e10;
            //will be used to erase nodes where all users come out
            long minIndex = 0;
            int i = -1, m = -1, t = -1;

            for(int z = 0; z < x; z++) {
                double curCost = beamSolution[z].cost;

                //if there are no users in the source cell (i), skip it
                if (problem.usersCell[beamSolution[z].i][beamSolution[z].m][beamSolution[z].t] == 0) continue;
                if(reminder != 0 && problem.n[beamSolution[z].m] > reminder) continue;
                if(j == beamSolution[z].i) continue;

                if(curCost < minCost || (curCost == minCost && beamSolution[z].m > m)) {
                    minCost = curCost;
                    minIndex = z;
                    i = beamSolution[z].i;
                    m = beamSolution[z].m;
                    t = beamSolution[z].t;
                }
            }

            if(problem.usersCell[i][m][t] * problem.n[m] > demand) {
                //not moving out all users
                reminder = demand % problem.n[m];
                int users_moving = demand/problem.n[m];
                problem.usersCell[i][m][t] -= users_moving;
                // move an user inside the destination
                //problem.usersCell[destination][m][t] += users_moving;
                solution[i][j][m][t] += users_moving;
                demand -= users_moving * problem.n[m];
            } else {
                //moving out all users
                // move an user inside the destination
                //problem.usersCell[destination][m][t] += problem.usersCell[i][m][t];
                solution[i][j][m][t] += problem.usersCell[i][m][t];
                // update number of users in source cell
                demand -= problem.n[m]*problem.usersCell[i][m][t];
                problem.usersCell[i][m][t] = 0;
                reminder = 0;

                beamSolution[minIndex].cost = 1000;
            }
        }
    }*/

    //calculate total number of tasks
    demand = 0;
    for(int i = 0; i < nCells; i++) {
        demand += problem.activities[i];
    }

    //initialize alleles
    GARealAlleleSet destSet;
    GARealAlleleSet sourceSet;
    for(int j = 0; j < nCells; j++) {
        if (problem.activities[j] != 0) {
            destSet.add(j);
        }
        else {
            sourceSet.add(j);
        }
    }
    GARealAlleleSet userTypeSet;
    for(int m = 0; m < nCustomerTypes; m++)
        userTypeSet.add(m);
    GARealAlleleSet timeSet;
    for(int t = 0; t < nTimeSteps; t++)
        timeSet.add(t);

    GARealAlleleSetArray alleles;
    alleles.add(sourceSet);
    alleles.add(destSet);
    alleles.add(userTypeSet);
    alleles.add(timeSet);

    GARealGenome genome(alleles, Objective);
    //genome.initializer(greedyInitializer);

    GASteadyStateGA ga(genome);
    //ga.nPopulations(5);
    ga.populationSize((unsigned int) 5*demand);
    //ga.pMigration((unsigned int) 0.1*ga.populationSize());
    ga.pReplacement(pRep != -1 ? pRep : 0.48);
    ga.scaling(GANoScaling()); // fitness value must be equal to the objective function
    ga.selector(GATournamentSelector()); //// pick always the best genome
    ga.pMutation(pMut != -1 ? pMut : 0.42);
    ga.pCrossover(pCross != -1 ? : 0.34);
    ga.nBestGenomes((unsigned int) demand); // find at least "demands" solutions

    ga.terminator(Heuristic::GATermination);
    ga.evolve();

    //const GAGenome & genome =  ga.statistics().bestPopulation().select();
    //ga.statistics().nBestGenomes(ga.statistics().bestPopulation().select(), ga.nBestGenomes());
    //ga.statistics().bestPopulation().order(GAPopulation::SortOrder::HIGH_IS_BEST);
    ga.statistics().bestPopulation().sort();

    // Update problem variables according to the found solution
    int i, j, m, t;
    for(unsigned int k = 0; k < ga.statistics().bestPopulation().size() && demand > 0; k++) {
        GARealGenome& g = (GARealGenome&) ga.statistics().bestIndividual(k);

        i = (int) g.gene(0);
        j = (int) g.gene(1);
        m = (int) g.gene(2);
        t = (int) g.gene(3);

        //std::cout << problem.costs[i][j][m][t]/problem.n[m] << std::endl;

        if(problem.usersCell[i][m][t] == 0) continue;
        if(problem.activities[j] == 0) continue;
        if(i == j) continue;


        int possibleTasks = problem.usersCell[i][m][t]*problem.n[m];


        if(problem.activities[j] >= possibleTasks) {
            solution[i][j][m][t] += problem.usersCell[i][m][t];
            problem.usersCell[i][m][t] = 0;
            demand -= possibleTasks;
            problem.activities[j] -= possibleTasks;
        } else {
            solution[i][j][m][t] += problem.activities[j] / problem.n[m];
            problem.usersCell[i][m][t] -= problem.activities[j] / problem.n[m];
            demand -= (problem.activities[j] / problem.n[m]) * problem.n[m];
            problem.activities[j] -= (problem.activities[j] / problem.n[m]) * problem.n[m];
        }

    }

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++) {
                    if(solution[i][j][m][t] > 0)
                        objFun += solution[i][j][m][t] * problem.costs[i][j][m][t];
                }

    stat.push_back(execTime);
    stat.push_back(objFun);

    hasSolution=true;

    return objFun;
}

float Heuristic::Objective(GAGenome & g) {
    GARealGenome& genome = (GARealGenome&) g;

    int i, j, m, t;
    i = (int) genome.gene(0);
    j = (int) genome.gene(1);
    m = (int) genome.gene(2);
    t = (int) genome.gene(3);

    int users = problem.usersCell[i][m][t];
    if(users == 0) return 0;
    // if there is nothing to do
    if(problem.activities[j] == 0) return 0;
    if(i == j) return 0;
    if(problem.activities[j] < problem.n[m]) return 0;

    return ((float)problem.n[m] / ((float)problem.costs[i][j][m][t])); // min cost? -> max 1/cost
}

GABoolean
Heuristic::GATermination(GAGeneticAlgorithm & ga) {
    execTime = (clock() - execTimeStart) / CLOCKS_PER_SEC;
    return(execTime < 4.8) ? gaFalse : gaTrue;
}

void Heuristic::greedyInitializer(GAGenome &g) {

    GARealGenome &genome=(GARealGenome &)g;
    cout << genome << endl;
    /*int z = 0;
    for (int i = 0; i < 30; i++)
        for (int j = 0; j < 30; j++)
            for (int m = 0; m < 3; m++)
                for (int t = 0; t < 1; t++, z++)
                    if(i != j || problem.activities[j] > 0);
                        // add gene to genoma*/

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


void Heuristic::upperHeuristic(int j, int m, int taskNumber, double timeLimit, PartialSolution* sol) {
    int x = 0;
    //i source cell
    for(int i = 0; i < nCells/* && x < taskNumber*/; i++) {
        for(int t = 0; t < nTimeSteps/* && x < taskNumber*/; t++, x++) {
            PartialSolution s;
            s.cost = problem.costs[i][j][m][t]/problem.n[m];
            s.i = i;
            s.m = m;
            s.t = t;

            sol[x] = s;
        }
    }
}