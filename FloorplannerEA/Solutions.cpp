#include "Solutions.h"

Solutions::Solutions(int anum_of_solutions,int anum_objectives, int anum_variables,int aL, int aW, int aH) {
    this->number_of_solutions = 2*anum_of_solutions;
    number_of_elements = anum_variables;
    number_of_objectives = anum_objectives;
    L = aL;
    W = aW;
    H = aH;
    ids = new int[number_of_elements];
    names = new string[number_of_elements];
    for(int i=0;i<number_of_elements;i++){
            names[i]="";
    }
    ps = new float[number_of_elements];
    dps = new float[number_of_elements];
    ls = new int[number_of_elements];
    ws = new int[number_of_elements];
    rotated = new bool*[number_of_solutions];
    rotated[0] = new bool [number_of_solutions*number_of_elements];
    for (int i = 1; i < number_of_solutions; ++i){
        rotated[i] = rotated[i-1] + number_of_elements;
    }
    pos = new int*[number_of_solutions];
    pos[0] = new int [number_of_solutions*number_of_elements];
    for (int i = 1; i < number_of_solutions; ++i){
        pos[i] = pos[i-1] + number_of_elements;
    }
    floorplan = new int* [number_of_solutions];
    floorplan[0] = new int [number_of_solutions*L*W*H];
    for (int i = 1; i < number_of_solutions; ++i){
        floorplan[i] = floorplan[i-1] + L*W*H;
    }
    surface = new int* [number_of_solutions];
    surface[0] = new int [number_of_solutions*L*W*H];
    for (int i = 1; i < number_of_solutions; ++i){
        surface[i] = surface[i-1] + L*W*H;
    }
    /*ls = new int* [number_of_solutions];
    ls[0] = new int[number_of_solutions*number_of_elements];
    for (int i = 1; i < number_of_solutions; ++i){
        ls[i] = ls[i-1] + number_of_elements;
    }
    ws = new int* [number_of_solutions];
    ws[0] = new int[number_of_solutions*number_of_elements];
    for (int i = 1; i < number_of_solutions; ++i){
        ws[i] = ws[i-1] + number_of_elements;
    }*/
    objectives = new float* [number_of_solutions];
    objectives[0] = new float[number_of_solutions*number_of_objectives];
    for (int i = 1; i < number_of_solutions; ++i){
        objectives[i] = objectives[i-1] + number_of_objectives;
    }
    crowdingDistances = new float[number_of_solutions];
    ranks = new int[number_of_solutions];
    rankValues = new int[number_of_solutions];
}

Solutions::~Solutions() {
    for(int i=0; i<number_of_elements;i++){
        names[i]="";
    }
    delete[] ids; ids = NULL;
    delete[] ls; ls = NULL;
    delete[] ws; ws = NULL;
    delete[] rotated; rotated = NULL;
    
    delete []dps; dps = NULL;
    delete[] pos; pos = NULL;
    delete[] floorplan;
    delete[] surface;
    delete [] objectives;
    delete ranks; ranks = NULL;
    delete rankValues; rankValues = NULL;
    delete crowdingDistances;
}


void Solutions::rankSolutions(){
    for(int i=0; i<number_of_solutions; i++) {
        ranks[i] = 0;
        rankValues[i] = number_of_solutions+1;
    }
    for(int i=0; i<number_of_solutions; i++) {
        for(int j=0; j<number_of_solutions; j++) {
            if(i!=j){
                if(compareSolutionDominance(i, j)>0) { // j dominates i
                    ranks[i]++;
                }
            }
        }
    }
    for(int i=0; i<number_of_solutions-1; i++) {
        int index = 0;
        while(ranks[i]>rankValues[index]){
            index++;
        }
        if(ranks[i] != rankValues[index]){
            insertRank(index,ranks[i]);
        }
    }
}

void Solutions::insertRank(int index, int rank){
    for(int i=number_of_solutions;i>index+1;i--){
        rankValues[i-1] = rankValues[i-2];
    }
    rankValues[index] = rank;
}

int Solutions::compareSolutionDominance(int indice_i, int indice_j){
    bool bigger = false;
    bool smaller = false;
    bool indiff = false;
    for(int i=0; !(indiff) && i<number_of_objectives; i++) {
        if(objectives[indice_i][i] > objectives[indice_j][i]) bigger  = true;
        if(objectives[indice_i][i] < objectives[indice_j][i]) smaller = true;
        indiff = (bigger && smaller);
    }
    if(smaller && !bigger) return -1;
    else if(bigger && !smaller) return 1;
    return 0;
}

int Solutions::GetIndexFromId(int indice_sol,int id){
    int index = pos[indice_sol][id-1];
    return index;
}

void Solutions::getCoordinates(int* CIx,int* CIy,int* CIz,int posI){
    int totalsize = this->L * this->W * this->H;
    int cellsPerLayer = totalsize / this->H;
    *CIz = (posI / cellsPerLayer);
    posI = posI - (cellsPerLayer * *CIz);
    *CIx = posI % this->L;
    posI = posI - *CIx;
    *CIy = posI / this->L;
}

void Solutions::replace(int newPos,int oldPos,int newRank){
    float* arrDouAux;
    arrDouAux = objectives[newPos];
    objectives[newPos] = objectives[oldPos];
    objectives[oldPos] = arrDouAux;
    bool* arrBoolAux;
    arrBoolAux = rotated[newPos];
    rotated[newPos] = rotated[oldPos];
    rotated[oldPos] = arrBoolAux;
    
    /*arrIntAux = ls[newPos];
    ls[newPos] = ls[oldPos];
    ls[oldPos] = arrIntAux;
    arrIntAux = ws[newPos];
    ws[newPos] = ws[oldPos];
    ws[oldPos] = arrIntAux;*/
    int* arrIntAux;
    arrIntAux = pos[newPos];
    pos[newPos] = pos[oldPos];
    pos[oldPos] = arrIntAux;
    arrIntAux = floorplan[newPos];
    floorplan[newPos] = floorplan[oldPos];
    floorplan[oldPos] = arrIntAux;
    arrIntAux = surface[newPos];
    surface[newPos] = surface[oldPos];
    surface[oldPos] = arrIntAux;
    float daux = crowdingDistances[newPos];
    crowdingDistances[newPos] = crowdingDistances[oldPos];
    crowdingDistances[oldPos] = daux;
    int iaux = ranks[newPos];
    //int iaux2 = ranks[oldPos];
    ranks[newPos] = ranks[oldPos];
    ranks[oldPos] = iaux;
    ranks[newPos] = newRank;
}

bool Solutions::feasiblePosition(int indSol,int id, int index_corner){
    int feasible=0;
    int l,w;
    if(!rotated[indSol][id-1]){
        l = ls[id-1];
        w = ws[id-1];
    }else{
        l = ws[id-1];
        w = ls[id-1];
    }
    int h = 1;
    int x,y,z;
    getCoordinates(&x,&y,&z,index_corner);
    if(x+l>L)feasible++;
    if(y+w>W)feasible++;
    if(z+h>H)feasible++;
    if(feasible==0){
        return true;
    }else{
        return false;
    }
}

void Solutions::solutionToXml(int indSol) {
    stringstream output;
    int cellSize=600;
    output  << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n" ;
    output  << "<Floorplan Version=\"5.0\" CellSize=\"" << cellSize
            << "\" Length=\"" << L << "\" Width=\"" << W
            << "\" NumLayers=\"" << H
            << "\" NumPowerProfiles=\"" << 1
            << "\">\n";
    output  << "\t<Blocks>\n";
    for(int i=0;i<number_of_elements;i++){
        int id=ids[i];
        int x,y,z,l,w;
        int pos= this->pos[indSol][id-1];
        getCoordinates(&x,&y,&z,pos);
        if(!rotated[indSol][id-1]){
            l = ls[id-1];
            w = ws[id-1];
        }else{
            l = ws[id-1];
            w = ls[id-1];
        }
        output  << "\t\t<Block id=\"" <<  id
                << "\" name=\"" << names[id-1]
                << "\" xMin=\"0\" x=\"" << x
                << "\" xMax=\"" << L - l
                << "\" yMin=\"0\" y=\"" << y
                << "\" yMax=\"" << W - w
                << "\" zMin=\"0\" z=\"" << z
                << "\" zMax=\"" << H - 1
                << "\" l=\"" << l
                << "\" w=\"" << w
                << "\" h=\"" << 1
                << "\" p=\"" << ps[id-1] << "\"/>\n";
                //<< "\" dp0=\"" << dps[id-1] << "\"/>\n";
    }
    output << "\t</Blocks>\n";
    output << "\t<Distances>\n";
    output << "\t</Distances>\n";
    output << "\t<Couplings>\n";
    output << "\t</Couplings>\n";
    output << "</Floorplan>\n";
    ofstream outFile;
    string name="";
    stringstream ssname;
    ssname << "results/" << indSol << ".xml";
    name = ssname.str();
    outFile.open(name.c_str());
    outFile << output.str();
}
