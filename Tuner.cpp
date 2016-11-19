//
// Created by mario on 11/11/16.
//

#include "heuristic.h"
#include <cmath>

float obj(GAGenome & g);
GABoolean termination(GAGeneticAlgorithm & ga);

int counter = 0;
float oldObj = -1;

#define N_INST 9
string INSTANCES[] = {
        /*"./input/Co_30_20_NT_0-txt", */
        "./input/Co_100_1_NT_0.txt",
        /*"./input/Co_100_1_NT_1.txt",
        "./input/Co_100_1_NT_2.txt",
        "./input/Co_100_1_T_0.txt"
        "./input/Co_100_1_T_1.txt"
        "./input/Co_100_1_T_2.txt",*/
        /*"./input/Co_300_20_NT_0.txt"*/
        //"./input/Co_30_1_NT_0.txt"
        /*"./input/Co_30_1_NT_1.txt",
        "./input/Co_30_1_NT_2.txt"*/
};
float OPTIMUM_SOLUTIONS[] = {
        /*719,*/
        4453,
        /*4530,
        4267,
        4276,
        4815,
        4636,*/
        /*7369*/
        //1041
        /*1756,
        2341*/
};

float obj(GAGenome &g) {
    GARealGenome& genome = (GARealGenome&) g;

    float pRep = genome.gene(1);
    float pMut = genome.gene(2);
    float pCross = genome.gene(3);

    Heuristic _heuristic(INSTANCES[counter / 10]);
    vector<double> stat;
    double obj = _heuristic.solveFast(stat, 5, pRep, pMut, pCross);

    cout << genome << endl;

    float opt = OPTIMUM_SOLUTIONS[counter / 10];
    if(opt > obj) {
        if(genome.geneticAlgorithm()->pMutation() < 0.9) genome.geneticAlgorithm()->pMutation(genome.geneticAlgorithm()->pMutation() + 0.1);
        return 0.0;
    }

    if(oldObj == -1) oldObj = obj;

    if(oldObj < obj) return 0.0;
    else oldObj = obj;

    counter++;

    double optimalityGap = (obj-opt)/opt;
    if(optimalityGap > 0.02) {
        if(genome.geneticAlgorithm()->pMutation() < 0.9) genome.geneticAlgorithm()->pMutation(genome.geneticAlgorithm()->pMutation() + 0.1);
        return 0.0;
    } else
        genome.geneticAlgorithm()->pMutation(0.01);

    cout << "Found: " << obj << " but optimum is: " << opt << " optimality gap is: " << optimalityGap << endl;
    cout << genome << endl;

    return 1 / abs(obj-opt);
}

GABool termination(GAGeneticAlgorithm &ga) {
    return counter >= N_INST*10 ? GABoolean::gaTrue : GABoolean::gaFalse ;
}

int main(int argc, char **argv) {
    // Genetic Algorithm
    GARealAlleleSet set;
    for(float i = 0.01; i < 1.0; i += 0.01)
        set.add(i);

    GARealAlleleSetArray alleles;
    alleles.add(set);
    alleles.add(set);
    alleles.add(set);

    GARealGenome genome(alleles, obj);

    GASteadyStateGA ga(genome);
    //ga.nPopulations(5);
    ga.populationSize((unsigned int) 20);
    //ga.pMigration((unsigned int) 0.1*ga.populationSize());
    ga.pReplacement(0.4);
    ga.scaling(GANoScaling()); // fitness value must be equal to the objective function
    ga.selector(GATournamentSelector()); //// pick always the best genome
    ga.pMutation(1.0);
    ga.pCrossover(0.5);
    ga.nBestGenomes((unsigned int) 5);

    //ga.nGenerations(10);
    ga.terminator(termination);

    ga.evolve();

    std::cout << INSTANCES[counter / 10] << std::endl;
    std::cout << ga.statistics().bestIndividual() << std::endl;

    return 0;
}
