#include <iostream>
#include <random>
#include <thread>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <sys/time.h>
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
    //initializing variables of the problem
    double objFun = 0;

    timeval tim;
    gettimeofday(&tim, NULL);
    execTimeStart = (tim.tv_sec*1000+tim.tv_usec/1000);

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;


    //definition and initialization of my own variables
    //objfun
    double *objfuns = new double[4];
    for(int i = 0; i < 4; i++) {
        objfuns[i] = 0;
    }

    //solution
    int ****solution1 = new int***[nCells];
    int ****solution2 = new int***[nCells];
    int ****solution3 = new int***[nCells];
    int ****solution4 = new int***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        solution1[i] = new int**[nCells];
        solution2[i] = new int**[nCells];
        solution3[i] = new int**[nCells];
        solution4[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            solution1[i][j] = new int*[nCustomerTypes];
            solution2[i][j] = new int*[nCustomerTypes];
            solution3[i][j] = new int*[nCustomerTypes];
            solution4[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                solution1[i][j][m] = new int[nTimeSteps];
                solution2[i][j][m] = new int[nTimeSteps];
                solution3[i][j][m] = new int[nTimeSteps];
                solution4[i][j][m] = new int[nTimeSteps];
            }
        }
    }
    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++) {
                    solution1[i][j][m][t] = 0;
                    solution2[i][j][m][t] = 0;
                    solution3[i][j][m][t] = 0;
                    solution4[i][j][m][t] = 0;
                }

    //usercell
    int ***usersCell1 = new int**[nCells];
    int ***usersCell2 = new int**[nCells];
    int ***usersCell3 = new int**[nCells];
    int ***usersCell4 = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        usersCell1[i] = new int*[nCustomerTypes];
        usersCell2[i] = new int*[nCustomerTypes];
        usersCell3[i] = new int*[nCustomerTypes];
        usersCell4[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            usersCell1[i][m] = new int[nTimeSteps];
            usersCell2[i][m] = new int[nTimeSteps];
            usersCell3[i][m] = new int[nTimeSteps];
            usersCell4[i][m] = new int[nTimeSteps];
        }
    }
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            for (int i = 0; i < nCells; i++) {
                usersCell1[i][m][t] = problem.usersCell[i][m][t];
                usersCell2[i][m][t] = problem.usersCell[i][m][t];
                usersCell3[i][m][t] = problem.usersCell[i][m][t];
                usersCell4[i][m][t] = problem.usersCell[i][m][t];
            }
        }
    }

    //activities
    int *activities1 = new int[nCells];
    int *activities2 = new int[nCells];
    int *activities3 = new int[nCells];
    int *activities4 = new int[nCells];
    for (int i = 0; i < nCells; i++) {
        activities1[i] = problem.activities[i];
        activities2[i] = problem.activities[i];
        activities3[i] = problem.activities[i];
        activities4[i] = problem.activities[i];
    }

    //demand
    int demand = 0;
    for(int i = 0; i < nCells; i++) {
        demand += problem.activities[i];
    }



    thread t1(greedy1, demand, activities1, usersCell1, solution1, objfuns);
    thread t2(greedy2, demand, activities2, usersCell2, solution2, objfuns);
    thread t3(greedy3, demand, activities3, usersCell3, solution3, objfuns);
    thread t4(greedy4, demand, activities4, usersCell4, solution4, objfuns);

    t1.join();
    t2.join();
    t3.join();
    t4.join();


    double min = objfuns[0];
    int index = 0;
    for(int x = 0; x < 4; x++) {
        if(objfuns[x] < min) {
            min = objfuns[x];
            index = x;
        }
    }


    switch (index) {
        case 0: {
            for (int i = 0; i < nCells; i++)
                for (int j = 0; j < nCells; j++)
                    for (int m = 0; m < nCustomerTypes; m++)
                        for (int t = 0; t < nTimeSteps; t++)
                            if (solution1[i][j][m][t] > 0) {
                                solution[i][j][m][t] = solution1[i][j][m][t];
                                problem.usersCell[i][m][t] -= solution[i][j][m][t];
                                problem.activities[j] = 0;
                                objFun += solution[i][j][m][t] * problem.costs[i][j][m][t];
                            }
        }
            break;
        case 1: {
            for (int i = 0; i < nCells; i++)
                for (int j = 0; j < nCells; j++)
                    for (int m = 0; m < nCustomerTypes; m++)
                        for (int t = 0; t < nTimeSteps; t++)
                            if (solution2[i][j][m][t] > 0) {
                                solution[i][j][m][t] = solution2[i][j][m][t];
                                problem.usersCell[i][m][t] -= solution[i][j][m][t];
                                problem.activities[j] = 0;
                                objFun += solution[i][j][m][t] * problem.costs[i][j][m][t];
                            }
        }
            break;
        case 2: {
            for (int i = 0; i < nCells; i++)
                for (int j = 0; j < nCells; j++)
                    for (int m = 0; m < nCustomerTypes; m++)
                        for (int t = 0; t < nTimeSteps; t++)
                            if (solution3[i][j][m][t] > 0) {
                                solution[i][j][m][t] = solution3[i][j][m][t];
                                problem.usersCell[i][m][t] -= solution[i][j][m][t];
                                problem.activities[j] = 0;
                                objFun += solution[i][j][m][t] * problem.costs[i][j][m][t];
                            }
        }
            break;
        case 3: {
            for (int i = 0; i < nCells; i++)
                for (int j = 0; j < nCells; j++)
                    for (int m = 0; m < nCustomerTypes; m++)
                        for (int t = 0; t < nTimeSteps; t++)
                            if (solution4[i][j][m][t] > 0) {
                                solution[i][j][m][t] = solution4[i][j][m][t];
                                problem.usersCell[i][m][t] -= solution[i][j][m][t];
                                problem.activities[j] = 0;
                                objFun += solution[i][j][m][t] * problem.costs[i][j][m][t];
                            }
        }
            break;
    }

    gettimeofday(&tim, NULL);
    execTime = (tim.tv_sec*1000+tim.tv_usec/1000) - execTimeStart;

    stat.push_back(execTime);
    stat.push_back(objFun);

    hasSolution = true;

    return objFun;
}

void Heuristic::greedy1(int demand1, int *activities1, int ***usersCell1, int ****solution1, double *obj) {
    //solution
    int it2 =-1, jt2 =-1, mt2 =-1, tt2=-1;
    while(demand1 > 0) {
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;

        for (int j = 0; j < nCells; j++) {
            if(activities1[j] == 0) continue;
            for (int i = 0; i < nCells; i++) {
                if(i == j || activities1[j] == 0) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    //if(activities[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(usersCell1[i][m][t]==0) continue;

                        if(minCost > ((double)problem.costs[i][j][m][t]*activities1[j]/(double)problem.n[m])
                           || (minCost == ((double)problem.costs[i][j][m][t]*activities1[j]/(double)problem.n[m]) && m > mt)) {
                            minCost = ((double)problem.costs[i][j][m][t]*activities1[j]/(double)problem.n[m]);
                            it = i; mt = m; tt = t; jt = j;
                        }
                    }
                }
            }
        }

        minCost = 1e10;
        if(activities1[jt] < problem.n[mt]) {
            for (int i = 0; i < nCells; i++) {
                if(i == jt) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    //if(activities[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(usersCell1[i][m][t]==0) continue;

                        if(minCost > ((double)problem.costs[i][jt][m][t]*activities1[jt])
                           || (minCost == problem.costs[i][jt][m][t]*activities1[jt] && m > mt)) {
                            minCost = ((double)problem.costs[i][jt][m][t]*activities1[jt]);
                            it = i; mt = m; tt = t;
                        }
                    }
                }
            }
        }

        // update solution1s
        int maxTasks = problem.n[mt] * usersCell1[it][mt][tt];

        if(maxTasks > activities1[jt]) {
            int movedUsers = activities1[jt] / problem.n[mt];
            if(it2 == it && mt2 == mt && jt2 == jt && tt2 == tt) {
                movedUsers = 1;
                demand1 -= activities1[jt];
                activities1[jt] = 0;
            }
            else {
                demand1 -= movedUsers * problem.n[mt];
                activities1[jt] -= movedUsers * problem.n[mt];
            }
            solution1[it][jt][mt][tt] += movedUsers;
            usersCell1[it][mt][tt] -= movedUsers;
            it2 =it; mt2 =mt; jt2=jt; tt2=tt;
            //cout << i << " " << j << " " << m << " " << t << endl;
        }
        else {
            solution1[it][jt][mt][tt] += usersCell1[it][mt][tt];
            activities1[jt] -= maxTasks;
            usersCell1[it][mt][tt] = 0;
            demand1 -= maxTasks;
        }
    }

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++) {
            if(i == j) continue;
            if(problem.activities[j] == 0) continue;
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if (solution1[i][j][m][t] > 0) obj[0] += solution1[i][j][m][t] * problem.costs[i][j][m][t];
        }
}

void Heuristic::greedy2(int demand2, int *activities2, int ***usersCell2, int ****solution2, double *obj) {
    //start greedy2
    int it2 =-1, jt2 =-1, mt2 =-1, tt2=-1;
    while(demand2 > 0) {
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;


        for (int j = 0; j < nCells; j++) {
            if(activities2[j] == 0) continue;
            for (int i = 0; i < nCells; i++) {
                if(i == j) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    //if(activities[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(usersCell2[i][m][t]==0) continue;

                        if(minCost > ((double)problem.costs[i][j][m][t]/(double)problem.n[m])
                           || (minCost == ((double)problem.costs[i][j][m][t]/(double)problem.n[m]) && m > mt)) {
                            minCost = ((double)problem.costs[i][j][m][t]/(double)problem.n[m]);
                            it = i; mt = m; tt = t; jt = j;
                        }
                    }
                }
            }
        }

        minCost = 1e10;
        if(activities2[jt] < problem.n[mt]) {
            for (int i = 0; i < nCells; i++) {
                if(i == jt) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    //if(activities[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(usersCell2[i][m][t]==0) continue;

                        if(minCost > ((double)problem.costs[i][jt][m][t])
                           || (minCost == problem.costs[i][jt][m][t] && m > mt)) {
                            minCost = ((double)problem.costs[i][jt][m][t]);
                            it = i; mt = m; tt = t;
                        }
                    }
                }
            }
        }

        // update solution1s
        int maxTasks = problem.n[mt] * usersCell2[it][mt][tt];

        if(maxTasks > activities2[jt]) {
            int movedUsers = activities2[jt] / problem.n[mt];

            if(it2 == it && mt2 == mt && jt2 == jt && tt2 == tt) {
                movedUsers = 1;
                demand2 -= activities2[jt];
                activities2[jt] = 0;
            }
            else {
                demand2 -= movedUsers * problem.n[mt];
                activities2[jt] -= movedUsers * problem.n[mt];
            }

            solution2[it][jt][mt][tt] += movedUsers;
            usersCell2[it][mt][tt] -= movedUsers;
            it2 =it; mt2 =mt; jt2=jt; tt2=tt;
            //cout << it << " " << jt << " " << mt << " " << tt << demand << " " << activities[jt] << endl;
        }
        else {
            solution2[it][jt][mt][tt] += usersCell2[it][mt][tt];
            activities2[jt] -= maxTasks;
            usersCell2[it][mt][tt] = 0;
            demand2 -= maxTasks;
        }
    }
    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++) {
            if(i == j) continue;
            if(problem.activities[j] == 0) continue;
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if (solution2[i][j][m][t] > 0) obj[1] += solution2[i][j][m][t] * problem.costs[i][j][m][t];
        }
}

void Heuristic::greedy3(int demand3, int *activities3, int ***usersCell3, int ****solution3, double *obj) {


    //start greedy 3
    while(demand3 > 0) {
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;

        for (int j = 0; j < nCells; j++) {
            if(activities3[j] == 0) continue;
            for (int i = 0; i < nCells; i++) {
                if(i == j || activities3[j] == 0) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    if(activities3[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(usersCell3[i][m][t]==0) continue;

                        if(minCost > ((double)problem.costs[i][j][m][t]*activities3[j]/(double)problem.n[m])
                           || (minCost == ((double)problem.costs[i][j][m][t]*activities3[j]/(double)problem.n[m]) && m > mt)) {
                            minCost = ((double)problem.costs[i][j][m][t]*activities3[j]/(double)problem.n[m]);
                            it = i; mt = m; tt = t; jt = j;
                        }
                    }
                }
            }
        }

        // update solution1s
        int maxTasks = problem.n[mt] * usersCell3[it][mt][tt];

        if(maxTasks > activities3[jt]) {
            int movedUsers = activities3[jt] / problem.n[mt];
            solution3[it][jt][mt][tt] += movedUsers;
            activities3[jt] -= movedUsers * problem.n[mt];
            usersCell3[it][mt][tt] -= movedUsers;
            demand3 -= movedUsers * problem.n[mt];
            //cout << i << " " << j << " " << m << " " << t << endl;
        }
        else {
            solution3[it][jt][mt][tt] += usersCell3[it][mt][tt];
            activities3[jt] -= maxTasks;
            usersCell3[it][mt][tt] = 0;
            demand3 -= maxTasks;
        }
    }

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++) {
            if(i == j) continue;
            if(problem.activities[j] == 0) continue;
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if (solution3[i][j][m][t] > 0) obj[2] += solution3[i][j][m][t] * problem.costs[i][j][m][t];
        }
}

void Heuristic::greedy4(int demand4, int *activities4, int ***usersCell4, int ****solution4, double *obj) {
    //start greedy 4
    while(demand4 > 0) {
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;

        for (int j = 0; j < nCells; j++) {
            if(activities4[j] == 0) continue;
            for (int i = 0; i < nCells; i++) {
                if(i == j || activities4[j] == 0) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    if(activities4[j] < problem.n[m]) continue;
                    for (int t = 0; t < nTimeSteps; t++) {
                        if(usersCell4[i][m][t]==0) continue;

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
        int maxTasks = problem.n[mt] * usersCell4[it][mt][tt];

        if(maxTasks > activities4[jt]) {
            int movedUsers = activities4[jt] / problem.n[mt];
            solution4[it][jt][mt][tt] += movedUsers;
            activities4[jt] -= movedUsers * problem.n[mt];
            usersCell4[it][mt][tt] -= movedUsers;
            demand4 -= movedUsers * problem.n[mt];
            //cout << i << " " << j << " " << m << " " << t << endl;
        }
        else {
            solution4[it][jt][mt][tt] += usersCell4[it][mt][tt];
            activities4[jt] -= maxTasks;
            usersCell4[it][mt][tt] = 0;
            demand4 -= maxTasks;
        }
    }
    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++) {
            if(i == j) continue;
            if(problem.activities[j] == 0) continue;
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if (solution4[i][j][m][t] > 0) obj[3] += solution4[i][j][m][t] * problem.costs[i][j][m][t];
        }
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