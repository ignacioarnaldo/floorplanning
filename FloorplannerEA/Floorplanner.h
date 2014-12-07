#include <iostream>
#include <vector>
#include <cmath>
#include "FloorplanConfiguration.h"
#include "Solutions.h"

using namespace std;

#ifndef FLOORPLANNER_H
#define	FLOORPLANNER_H

class Floorplanner {

    public:
        Floorplanner(){};
        ~Floorplanner();
        Floorplanner(string xmlFilePath,int aMaxPopulation, int aMaxGeneration,int anumObjectives);
        void initialize();
        void newRandomSetOfSolutions();
        void evaluate(int indice_i);
        void evaluateChildren();
        void computeTempAndWire(int indice_i,int length,float res[]);
        void computeTempAndWireChildren(float* objTemp, float* objWire);
	void computeWireTempKernelCPU(int maxPop,int numberOfVars,float* objT, float* objW,int* ids,float* dps,int* ls,
				      int* ws,bool** rotated,int** pos,bool** couplings,int H,int L,int W);
        int checkfeasibility(int indice_i);
        void execute();
        void step();
        void clearOffspring();
        int executeBinaryTournamentNSGAII();
        void SinglePointCrossover(int ind_parent1,int ind_parent2,int ind_child1,int ind_child2);
        void executeMutation(int Solind);
        void reduce(int maxSize);
        void saveNonDominatedFront();
        void printSolution(int indSol);
        void printFloorplan(int indSol);

        Solutions* population;
        int numberOfVariables, numberOfObjectives, maxPopulation, 
            maxGenerations,currentGeneration;
        float mutation_probability;
        int maxHeight, maxLength, maxWidth, cellSize;
        bool** couplings;
        FloorplanConfiguration cfg;
        float* lowerBound;
        float* upperBound;
        string name;
        float infinity;
        float bestObjs[2];
        vector<ComponentVariable> components;
};

#endif	/* FLOORPLANNER_H */
