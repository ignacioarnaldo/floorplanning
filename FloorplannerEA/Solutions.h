#include <fstream>
#include<sstream>

using namespace std;


#ifndef SOLUTIONS_H
#define	SOLUTIONS_H

class Solutions {
public:
    Solutions(){}
    Solutions(int anum_of_solutions,int anum_objectives, int anum_variables,int aL, int aW, int aH);
    ~Solutions();
    void rankSolutions();
    int GetIndexFromId(int indice_sol,int id);
    void getCoordinates(int* CIx,int* CIy,int* CIz,int posI);
    void insertRank(int index, int rank);
    int compareSolutionDominance(int indice_i, int indice_j);
    bool feasiblePosition(int indSol,int indComp, int index_corner);
    void replace(int newPos,int oldPos,int newRank);
    void solutionToXml(int indSol);
    int number_of_solutions;
    int number_of_elements, number_of_objectives, L, W, H;
    int* ids;
    string* names;
    float* ps;
    float* dps; // Power density
    int** floorplan;
    int** surface;
    int* ls;
    int* ws;
    bool** rotated;
    float** objectives;
    int* ranks;
    int* rankValues;
    int** pos;
    float* crowdingDistances;
    
private:
};

#endif	/* SOLUTIONS_H */
