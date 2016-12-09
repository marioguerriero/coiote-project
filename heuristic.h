#ifndef COIOTE_HEURISTIC_HEURISTIC_H
#define COIOTE_HEURISTIC_HEURISTIC_H

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <list>
#include <map>

#include <ga/ga.h>
#include <ga/GASimpleGA.h>	// we're going to use the simple GA
#include <ga/GAGenome.h>
#include <ga/GAStatistics.h>
#include <ga/GARealGenome.h>
#include <ga/std_stream.h>

using namespace std;

struct Data {
    /**
     * the costs ijkt
     */
    double**** costs;

    /**
     * number of activity done by one user of type k
     */
    int* n;

    /**
     * Activities to be done in node i during the given time horizon
     */
    int* activities;

    /**
     * Number of users of type m in i during time period t (ikt)
     */
    int*** usersCell;
};


struct PartialSolution {
    double cost;
    int i;
    int m;
    int t;
};


enum eFeasibleState {
    FEASIBLE,
    NOT_FEASIBLE_DEMAND,
    NOT_FEASIBLE_USERS
};

class Heuristic{
private:
    /**
     * Number of periods
     */
    static int nTimeSteps;

    /**
     * Number of customer types
     */
    static int nCustomerTypes;

    /**
     * Problem structure for parameters
     */
    static Data problem;

    /**
     * Flag equals to true if the problems has a solution
     */
    bool hasSolution;

    /**
     * Variables of the problem (X in the model)
     */
    static int**** solution;

    /**
     * How many tasks need to be done
     */
    static int demand;

    /**
     * Execution time
     */
    static double execTime;

    static double execTimeStart;

    static void greedy1(int demand1, int *activities1, int ***usersCell1, int ****solution1, double *obj);
    static void greedy2(int demand2, int *activities2, int ***usersCell2, int ****solution2, double *obj);
    static void greedy3(int demand3, int *activities3, int ***usersCell3, int ****solution3, double *obj);
    static void greedy4(int demand4, int *activities4, int ***usersCell4, int ****solution4, double *obj);

    static vector<int*> getUsersCombination(int demand);

public:
    /**
     * Default constructor
     * @return Heuristic object
     */
    Heuristic(){};

    /**
     * Constructor from external file
     * @param path path of the external file cotaining the instance of the problem
     * @return
     */
    Heuristic(string path);

    /**
     * Number of cells
     */
    static int nCells;

    /**
     * Function to CHANGE!!! This function only makes a very bad solution for the problem
     * @param stat Array of statistics. In position 0 put the objVal, in position 1 the computational time
     * @param timeLimit Time limit for computation
     * @param verbose
     */
    double solveFast(vector<double>& stat, int timeLimit = - 1, float pRep = -1, float pMut = -1, float pCross = -1);

    /**
     * Puts KPIs in the statistics' array. Call this only if problem has a solution
     * @param stat Array of statistics
     */
    void getStatSolution(vector<double>& stat);

    /**
     * Write KPIs on a file
     * @param path path of the file
     * @param nameInstance name of the instance
     * @param stat array of statistics
     */
    void writeKPI(string path, string nameInstance, vector<double> stat);

    /**
     * Write the detailed solution on a file
     * @param path path of the file
     */
    void writeSolution(string path);

    /**
     * Check the feasibility of the problem
     * @param path path of the solution file
     * @return a state of the check (i.e. FEASIBLE if the solution is feasible)
     */
    eFeasibleState isFeasible(string path);

    static float Objective(GAGenome &g);

    static GABoolean GATermination(GAGeneticAlgorithm &ga);
};

#endif //COIOTE_HEURISTIC_HEURISTIC_H
