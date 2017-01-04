#include <iostream>
#include <random>
#include <thread>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <sys/time.h>
#include <unordered_set>
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
struct timespec Heuristic::start;

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


double Heuristic::solveFast(vector<double>& stat) {
    //initializing variables of the problem
    double objFun = 0;

    struct timespec finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    //definition and initialization of my own variables
    //objfun
    double *objfuns = new double[4];
    for(int i = 0; i < 4; i++) {
        objfuns[i] = 0;
    }

    int ****solution1 = new int***[nCells];
    int ****solution2 = new int***[nCells];
    int ****solution3 = new int***[nCells];
    int ****solution4 = new int***[nCells];

    int ***usersCell1 = new int**[nCells];
    int ***usersCell2 = new int**[nCells];
    int ***usersCell3 = new int**[nCells];
    int ***usersCell4 = new int**[nCells];

    int *activities1 = new int[nCells];
    int *activities2 = new int[nCells];
    int *activities3 = new int[nCells];
    int *activities4 = new int[nCells];

    for (int i = 0; i < this->nCells; i++) {
        solution1[i] = new int**[nCells];
        solution2[i] = new int**[nCells];
        solution3[i] = new int**[nCells];
        solution4[i] = new int**[nCells];
        usersCell1[i] = new int*[nCustomerTypes];
        usersCell2[i] = new int*[nCustomerTypes];
        usersCell3[i] = new int*[nCustomerTypes];
        usersCell4[i] = new int*[nCustomerTypes];
        activities1[i] = problem.activities[i];
        activities2[i] = problem.activities[i];
        activities3[i] = problem.activities[i];
        activities4[i] = problem.activities[i];
        demand += problem.activities[i];
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
                usersCell1[i][m] = new int[nTimeSteps];
                usersCell2[i][m] = new int[nTimeSteps];
                usersCell3[i][m] = new int[nTimeSteps];
                usersCell4[i][m] = new int[nTimeSteps];
                for (int t = 0; t < this->nTimeSteps; t++) {
                    solution1[i][j][m][t] = 0;
                    solution2[i][j][m][t] = 0;
                    solution3[i][j][m][t] = 0;
                    solution4[i][j][m][t] = 0;
                    solution[i][j][m][t] = 0;
                    usersCell1[i][m][t] = problem.usersCell[i][m][t];
                    usersCell2[i][m][t] = problem.usersCell[i][m][t];
                    usersCell3[i][m][t] = problem.usersCell[i][m][t];
                    usersCell4[i][m][t] = problem.usersCell[i][m][t];
                }
            }
        }
    }

    //best in order: 3, 1, 2
    //1 hard
    //4 best greedy
#ifdef TWO_THREADS
    thread t1(combined_search, demand, activities1, usersCell1, solution1, objfuns);
    thread t4(greedy3, demand, activities4, usersCell4, solution4, objfuns);
#elif THREE_THREADS
    thread t1(combined_search, demand, activities1, usersCell1, solution1, objfuns);
    thread t2(greedy1, demand, activities2, usersCell2, solution2, objfuns);
    thread t4(greedy3, demand, activities4, usersCell4, solution4, objfuns);
#else
    thread t1(combined_search, demand, activities1, usersCell1, solution1, objfuns);
    thread t2(greedy1, demand, activities2, usersCell2, solution2, objfuns);
    thread t3(greedy2, demand, activities3, usersCell3, solution3, objfuns);
    thread t4(greedy3, demand, activities4, usersCell4, solution4, objfuns);
#endif

#ifdef TWO_THREADS
    t1.join();
    t4.join();
#elif THREE_THREADS
    t1.join();
    t2.join();
    t4.join();
#else
    t1.join();
    t2.join();
    t3.join();
    t4.join();
#endif

    double min = 0;
    int index = -1;
    for(int x = 0; x < 4; x++) {
        if(!objfuns[x]) continue;
        if(min == 0 || objfuns[x] < min) {
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

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    stat.push_back(elapsed);
    stat.push_back(objFun);

    if(objFun != 0)
        hasSolution = true;

    return objFun;
}


/*
 * goal of this method is to find a feasible solution with the best outcome we can achive.
 * with the following we can solve almost all the hard instances too.
 *
 * this is divided into 3 main steps:
 * 1) managing "sure moves"
 * 2) find a feasible solution
 * 3) swap until there is time
 *
 * the main idea of the algorithm is to solve the problem by using "combinations":
 * a combination is a set of users which solve in the most precise way the given (cell) demand.
 * for example: if we have a cell with 6 tasks to do, and three types of users 3-6-10. we can solve
 * that problem by using one 6-tasks user or two 3-tasks users.
 *
 * step (1)
 * "sure users" are those which has to be moved. for example certain cells has a reduced number of tasks
 * and as a result only one possible combination of users; moreover, others have a set of combinations where
 * a particular type of user necessarily has to be moved.
 * so in this first step we search for those users and then for the less cost moves
 *
 * step (2)
 * generation of new combinations (since sure moves have modified the initial problem) of all cells
 * search for the min cost one, looking at the available users
 * apply the move
 *
 * step (2)
 * swap until 5 sec
 */
void Heuristic::combined_search(int demand1, int *activities1, int ***usersCell1, int ****solution1, double *obj) {
    /*
     * step (1)
     */
    vector<vector<int *>> combinations;
    for(int j = 0; j < nCells; j++) {
        int tempdemand = activities1[j];
        vector<int *> combs = getUsersCombination(tempdemand);
        while (tempdemand != 0 && combs.size() == 0)
            combs = getUsersCombination(++tempdemand);
        combinations.push_back(combs);
    }
    //tabulist used to skip already done moves
    std::unordered_set<int> tabulist;
    //external loop: search and perform min cost move
    for(int end = 0; end < nCells; end++) {
        int *besttype; //pointer to the best found
        vector<vector<int>> bestsources; //will contains the min sources (i, t) found

        //start searching for the min cost solution for the given j
        int minj = -1;
        double minfound = 1e10;
        for(int j = 0; j < nCells; j++) {
            if (activities1[j] == 0 || tabulist.count(j)) continue;

            //calc minimum users moved for every type (m)
            //second temp array used to memorize current moves, the first one will be broken
            int *type = new int[nCustomerTypes];
            int *typecopy = new int[nCustomerTypes];
            for (int m = 0; m < nCustomerTypes; m++) {
                type[m] = combinations[j][0][m];
                typecopy[m] = combinations[j][0][m];
            }
            for (int x = 0; x < combinations[j].size(); x++) {
                for (int m = 0; m < nCustomerTypes; m++) {
                    if (type[m] > combinations[j][x][m]) {
                        type[m] = combinations[j][x][m];
                        typecopy[m] = combinations[j][x][m];
                    }
                }
            }

            vector<vector<int>> bestsourcestemp;
            //calc cost of the move for the given j cell
            double movecost = 0;
            for (int m = 0; m < nCustomerTypes; m++) {
                IS_TIME_OVER
                double mincost = 1e10;
                int it, tt;
                vector<int> temp;
                while (type[m] != 0) {
                    for (int i = 0; i < nCells; i++) {
                        if (i == j) continue;
                        for (int t = 0; t < nTimeSteps; t++) {
                            if (usersCell1[i][m][t] == 0) continue;

                            if (mincost > problem.costs[i][j][m][t]) {
                                mincost = problem.costs[i][j][m][t];
                                it = i;
                                tt = t;
                            }
                        }
                    }
                    temp.push_back(it);
                    temp.push_back(tt);
                    //here mincost has been found
                    //let's calc movecost
                    if (usersCell1[it][m][tt] >= type[m]) {
                        usersCell1[it][m][tt] -= type[m];
                        movecost += mincost * type[m];
                        type[m] = 0;
                    } else {
                        //if here, another cycle has to be done
                        movecost += mincost * usersCell1[it][m][tt];
                        type[m] -= usersCell1[it][m][tt];
                        usersCell1[it][m][tt] = 0;
                        mincost = 1e10;
                    }
                }
                bestsourcestemp.push_back(temp);
            }
            //every j has to be evaluated with the same users in order to choose the minimum one
            //reload users
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    for (int i = 0; i < nCells; i++)
                        usersCell1[i][m][t] = problem.usersCell[i][m][t];

            if(movecost < minfound) {
                minj = j;
                minfound = movecost;
                besttype = typecopy;
                bestsources = bestsourcestemp;
            }
            IS_TIME_OVER
        }
        if(minj == -1) break;

        //here i, j and t of the best solutions (min combination cost) have been found
        tabulist.insert(minj);
        //apply minimum move
        for(int m = 0; m < nCustomerTypes; m++) {
            if(besttype[m] == 0) continue;

            vector<int> its = bestsources[m];
            int it, tt;
            for(int x = 0; x < its.size(); x++) {
                it = its[x];
                tt = its[++x];
                if(usersCell1[it][m][tt] >= besttype[m]) {
                    demand1 -= besttype[m] * problem.n[m];
                    activities1[minj] -= besttype[m] * problem.n[m];
                    solution1[it][minj][m][tt] += besttype[m];
                    usersCell1[it][m][tt] -= besttype[m];
                    problem.usersCell[it][m][tt] -= besttype[m];
                    besttype[m] = 0;
                }
                else {
                    solution1[it][minj][m][tt] += usersCell1[it][m][tt];
                    activities1[minj] -= usersCell1[it][m][tt]*problem.n[m];
                    demand1 -= usersCell1[it][m][tt]*problem.n[m];
                    besttype[m] -= usersCell1[it][m][tt];
                    usersCell1[it][m][tt] = 0;
                    problem.usersCell[it][m][tt] = 0;
                }

                if(activities1[minj] < 0) {
                    demand1 += -activities1[minj];
                    activities1[minj] = 0;
                }
            }
        }
        IS_TIME_OVER
    }

    /*
     * step (2)
     */
    combinations.clear();
    for(int j = 0; j < nCells; j++) {
        int tempdemand = activities1[j];
        vector<int *> combs = getUsersCombination(tempdemand);
        while (tempdemand != 0 && combs.size() == 0)
            combs = getUsersCombination(++tempdemand);
        combinations.push_back(combs);
    }

    while(demand1 > 0) {
        //will contain the chosen combination
        int *usertypecombin;
        int *availableusers;
        //calc the available users every time, because cells have to be evaluated equally
        availableusers = new int[nCustomerTypes];
        for (int m = 0; m < nCustomerTypes; m++) {
            availableusers[m] = 0;
            for (int i = 0; i < nCells; i++)
                for (int t = 0; t < nTimeSteps; t++)
                    availableusers[m] += problem.usersCell[i][m][t];
        }

        vector<vector<int>> bestsources; //will contains the min sources (i, t) found

        //start search of the best combinations above all cells
        double min = 1e10;
        int minj;
        for (int j = 0; j < nCells; j++) {
            if(activities1[j] == 0) continue;

            //find the best combination in combinations vector
            int *bestcomb;
            double mincombinationcost = 1e10;
            vector<vector<int>> bestsourcestemp2;
            for (int x = 0; x < combinations[j].size(); x++) {
                vector<vector<int>> bestsourcestemp;
                int *currComb = combinations[j][x];
                double combinationcost = 0;

                //temp array is a copy of currcomb, and is used to keep track of moved users in order to select the proper combination
                int *temp = new int[nCustomerTypes];
                for (int i = 0; i < nCustomerTypes; i++)
                    temp[i] = currComb[i];

                //calc cost of the current combination (temp contains it)
                bool fail = false;
                for (int m = 0; m < nCustomerTypes; m++) {
                    if (availableusers[m] < temp[m]) {
                        fail = true;
                        break;
                    }

                    //additional "while cycle" used to choose the right combination considering also users in the cells
                    double mcost = 1e10;
                    int it = -1, tt = -1;
                    vector<int> tempit;
                    while (temp[m] != 0) {
                        IS_TIME_OVER
                        for (int i = 0; i < nCells; i++) {
                            if (i == j) continue;
                            for (int t = 0; t < nTimeSteps; t++) {
                                if (usersCell1[i][m][t] == 0) continue;
                                if (problem.costs[i][j][m][t] < mcost) {
                                    mcost = problem.costs[i][j][m][t];
                                    it = i;
                                    tt = t;
                                }
                            }
                        }
                        tempit.push_back(it);
                        tempit.push_back(tt);
                        if (usersCell1[it][m][tt] > temp[m]) {
                            usersCell1[it][m][tt] -= temp[m];
                            combinationcost += mcost * temp[m];
                            temp[m] = 0;
                        } else {
                            combinationcost += mcost * usersCell1[it][m][tt];
                            temp[m] -= usersCell1[it][m][tt];
                            usersCell1[it][m][tt] = 0;
                            mcost = 1e10;
                        }
                    }
                    bestsourcestemp.push_back(tempit);
                }


                //reload users
                for (int m = 0; m < nCustomerTypes; m++)
                    for (int t = 0; t < nTimeSteps; t++)
                        for (int i = 0; i < nCells; i++)
                            usersCell1[i][m][t] = problem.usersCell[i][m][t];
                if (!fail && combinationcost < mincombinationcost) {
                    mincombinationcost = combinationcost;
                    bestcomb = currComb;
                    bestsourcestemp2 = bestsourcestemp;
                }
            }

            if(mincombinationcost/activities1[j] < min) {
                minj = j;
                min = mincombinationcost/activities1[j];
                usertypecombin = bestcomb;
                bestsources = bestsourcestemp2;
            }
            IS_TIME_OVER
        }

        //min cost solution has been found, apply moves
        for(int m = 0; m < nCustomerTypes; m++) {
            if(usertypecombin[m] == 0) continue;

            vector<int> its = bestsources[m];

            int it, tt;
            for(int x = 0; x < its.size(); x++) {
                it = its[x];
                tt = its[++x];

                if(problem.usersCell[it][m][tt] >= usertypecombin[m]) {
                    demand1 -= usertypecombin[m] * problem.n[m];
                    activities1[minj] -= usertypecombin[m] * problem.n[m];
                    solution1[it][minj][m][tt] += usertypecombin[m];
                    problem.usersCell[it][m][tt] -= usertypecombin[m];
                    usersCell1[it][m][tt] -= usertypecombin[m];
                    availableusers[m] -= usertypecombin[m];
                    usertypecombin[m] = 0;
                }
                else {
                    solution1[it][minj][m][tt] += problem.usersCell[it][m][tt];
                    activities1[minj] -= problem.usersCell[it][m][tt]*problem.n[m];
                    demand1 -= problem.usersCell[it][m][tt]*problem.n[m];
                    usertypecombin[m] -= problem.usersCell[it][m][tt];
                    availableusers[m] -= problem.usersCell[it][m][tt];
                    problem.usersCell[it][m][tt] = 0;
                    usersCell1[it][m][tt] = 0;
                }

                if(activities1[minj] < 0) {
                    demand1 += -activities1[minj];
                    activities1[minj] = 0;
                }
            }
        }
        IS_TIME_OVER
    }

    /*
     * step (3)
     */
    swapsolutions(solution1, usersCell1);


    for (int j = 0; j < nCells; j++) {
        if(problem.activities[j] == 0) continue;
        for (int i = 0; i < nCells; i++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if (solution1[i][j][m][t] > 0) obj[0] += solution1[i][j][m][t] * problem.costs[i][j][m][t];
    }
}


void Heuristic::greedy1(int demand2, int *activities2, int ***usersCell2, int ****solution2, double *obj) {
    //start greedy2
    int it2 =-1, jt2 =-1, mt2 =-1, tt2=-1;
    while(demand2 > 0) {
        IS_TIME_OVER
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
        } else {
            solution2[it][jt][mt][tt] += usersCell2[it][mt][tt];
            activities2[jt] -= maxTasks;
            usersCell2[it][mt][tt] = 0;
            demand2 -= maxTasks;
        }
    }

    swapsolutions(solution2, usersCell2);

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if(solution2[i][j][m][t] > 0) obj[1] += solution2[i][j][m][t] * problem.costs[i][j][m][t];
}

void Heuristic::greedy2(int demand3, int *activities3, int ***usersCell3, int ****solution3, double *obj) {
    //start greedy 3
    while(demand3 > 0) {
        IS_TIME_OVER
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;

        for (int i = 0; i < nCells; i++) {
            for (int j = 0; j < nCells; j++) {
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

    swapsolutions(solution3, usersCell3);

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if(solution3[i][j][m][t] > 0) obj[2] += solution3[i][j][m][t] * problem.costs[i][j][m][t];
}

void Heuristic::greedy3(int demand4, int *activities4, int ***usersCell4, int ****solution4, double *obj) {
    //start greedy 4
    while(demand4 > 0) {
        IS_TIME_OVER
        double minCost = 1e10;
        int it = 0, jt = 0, mt=0, tt=0;

        for (int i = 0; i < nCells; i++) {
            for (int j = 0; j < nCells; j++) {
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

    swapsolutions(solution4, usersCell4);


    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    if(solution4[i][j][m][t] > 0) obj[3] += solution4[i][j][m][t] * problem.costs[i][j][m][t];
}



/*
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
            IS_TIME_OVER
            if (problem.activities[j] == 0) continue;
            for (int i = 0; i < nCells; i++) {
                IS_TIME_OVER
                if (i == j) continue;
                for (int m = 0; m < nCustomerTypes; m++) {
                    IS_TIME_OVER
                    for (int t = 0; t < nTimeSteps; t++) {
                        IS_TIME_OVER
                        if (solutions[i][j][m][t] == 0) continue;

                        // i -> j cost
                        int ijcost = (int)problem.costs[i][j][m][t];

                        for(int m2 = 0; m2 < nCustomerTypes; m2++) {
                            IS_TIME_OVER
                            //those are to check if it is possible to swap two different type users
                            int quotient = problem.n[m2] / problem.n[m];
                            int residual = problem.n[m2] % problem.n[m];
                            //if users are of the same type (or a combinable one) residual is 0
                            if (residual != 0) continue;

                            if (quotient == 1) {
                                //quotient == 1 -> m2 == m: swap one to one.
                                //restart search of the max difference
                                for (int j2 = 0; j2 < nCells; j2++) {
                                    IS_TIME_OVER
                                    if (problem.activities[j2] == 0) continue;
                                    // i -> j2 cost
                                    int ij2cost = (int)problem.costs[i][j2][m][t];
                                    for (int i2 = 0; i2 < nCells; i2++) {
                                        IS_TIME_OVER
                                        for (int t2 = 0; t2 < nTimeSteps; t2++) {
                                            IS_TIME_OVER
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
                                    IS_TIME_OVER
                                    if (problem.activities[j2] == 0) continue;
                                    int ij2cost = (int)problem.costs[i][j2][m][t];
                                    for (int i2 = 0; i2 < nCells; i2++) {
                                        IS_TIME_OVER
                                        for (int t2 = 0; t2 < nTimeSteps; t2++) {
                                            IS_TIME_OVER
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
        //cout << "Sono qui" << endl;
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