#include "Floorplanner.h"

Floorplanner::Floorplanner(string xmlFilePath,int aMaxPopulation, int aMaxGenerations,int aNumObjectives){
    bestObjs[0]=FLT_MAX ;
    bestObjs[1]=FLT_MAX ;
    infinity=FLT_MAX ;
    maxPopulation = aMaxPopulation;
    maxGenerations = aMaxGenerations;
    numberOfObjectives = aNumObjectives;

    cfg = FloorplanConfiguration(xmlFilePath);
    cellSize = cfg.cellSize;
    maxLength = cfg.maxLength;
    maxWidth = cfg.maxWidth;
    maxHeight = cfg.maxHeight;
    components= vector<ComponentVariable>(numberOfVariables);

    vector<ComponentVariable> values;
    map<int,ComponentVariable> acomponents=cfg.components;
    this->numberOfVariables=acomponents.size();
    components= vector<ComponentVariable>(numberOfVariables);
    for (int i=0;i<numberOfVariables;i++){
        ComponentVariable cvAux = acomponents.at(i+1);
        components.at(i) = cvAux;
    }
    lowerBound=new float[numberOfVariables];
    for(int i=0;i<numberOfVariables;i++){
        this->lowerBound[i]=0.0;
    }
    upperBound=new float[numberOfVariables];
    for(int i=0;i<numberOfVariables;i++){
        this->upperBound[i]=0.0;
    }
    couplings = new bool* [numberOfVariables];
    couplings[0] = new bool[numberOfVariables*numberOfVariables];
    for (int i = 1; i < numberOfVariables; ++i){
        couplings[i] = couplings[i-1] + numberOfVariables;
    }
    for(int i=0;i<numberOfVariables;i++){
        for(int j=0;j<numberOfVariables;j++){
            int iaux = cfg.couplings[i][j];
            couplings[i][j] = iaux;
        }
    }
}

Floorplanner::~Floorplanner() {
    delete[] lowerBound;
    lowerBound=NULL;
    delete[] upperBound;
    upperBound=NULL;
    delete[] couplings;
    couplings = NULL;
    delete population;
    population = NULL;
}

void Floorplanner::initialize() {
    this->mutation_probability = 1.0 / numberOfVariables;
    this->population=new Solutions(this->maxPopulation,numberOfObjectives,
            numberOfVariables,maxLength,maxWidth,maxHeight);
    newRandomSetOfSolutions();
    int tam=population->number_of_solutions;
    for(int i=0; i<tam/2;i++){
        evaluate(i);
    }
    currentGeneration = 0;
}

void Floorplanner::newRandomSetOfSolutions() {
    int num_solutions = this->maxPopulation;
    int totalCells= maxHeight * maxWidth * maxLength;
    int index_corners=0;
    for (unsigned int i=0;i<components.size();i++){
        int idAux = components.at(i).id;
        population->ids[i]   = idAux;
        population->names[idAux-1] = components.at(i).name;
        //population->types[idAux-1] = components.at(i).type;
        //population->hs[j]    = cpv->getH();
        population->ps[idAux-1]    = components.at(i).p;
        population->dps[idAux-1]   = components.at(i).dp;
        population->ls[idAux-1] = components.at(i).l;
        population->ws[idAux-1] = components.at(i).w;
    }

    for (int s = 0; s < num_solutions; s++){
        for(int i=0;i<maxHeight*maxWidth*maxLength;i++){
            population->floorplan[s][i]=0;
        }
    }
    
    for (int s = 0; s < num_solutions; s++){
        index_corners=rand()%totalCells;
        for(int i=0;i<numberOfVariables;i++){
            int id=components.at(i).id;
            population->rotated[s][id-1] = false;
            while((population->floorplan[s][index_corners] !=0) || (!population->feasiblePosition(s,id,index_corners))){
                index_corners = rand() % totalCells;
            }
            population->floorplan[s][index_corners] = id;
            population->pos[s][id-1] = index_corners;
        }   
    }
}

void Floorplanner::evaluate(int indice_i) {
    float feasible = 0;
    feasible = checkfeasibility(indice_i);
    float objs[2];
    computeTempAndWire(indice_i,population->number_of_elements,objs);
    population->objectives[indice_i][0] = feasible;
    population->objectives[indice_i][1] = (feasible + 1) * objs[0];
    population->objectives[indice_i][2] = (feasible + 1) * objs[1];
    if (feasible == 0 && objs[0]<bestObjs[0] && objs[1]<bestObjs[1]){
        cout << "One feasible solution found: (" << feasible << "," << objs[0] << "," << objs[1] << ")" << endl;
        bestObjs[0] = objs[0];
        bestObjs[1] = objs[1];
    }
}

void Floorplanner::computeTempAndWire(int indice_sol,int length,float res[]){
    float tempObj = 0;
    float wireObj = 0;
    int CIz,CIy,CIx;
    int CJz,CJy,CJx;
    int CIid,posI,CIl,CIw,CIh;
    float CIdp;
    int CJid,posJ,CJl,CJw,CJh;
    float CJdp;
    for (int i = 0; i < length; ++i) {
        CIid = population->ids[i];
        posI = population->GetIndexFromId(indice_sol,CIid);
        population->getCoordinates(&CIx,&CIy,&CIz,posI);
        //CIl = population->ls[indice_sol][CIid-1];
        if(!population->rotated[indice_sol][CIid-1]){
            CIl = population->ls[CIid-1];
            CIw = population->ws[CIid-1];
        }else{
            CIl = population->ws[CIid-1];
            CIw = population->ls[CIid-1];
        }
        //int CIh = solutions->hs[indice_sol][i];
        CIh = 1;
        CIdp = population->dps[CIid-1];
        for (int j = 0; j < length; ++j) {
            if (i < j) {
                CJid = population->ids[j];
                posJ = population->GetIndexFromId(indice_sol,CJid);
                population->getCoordinates(&CJx,&CJy,&CJz,posJ);
                if(!population->rotated[indice_sol][CJid-1]){
                    CJl = population->ls[CJid-1];
                    CJw = population->ws[CJid-1];
                }else{
                    CJl = population->ws[CJid-1];
                    CJw = population->ls[CJid-1];
                }
                //int CJh = solution->hs[j];
                CJh = 1;
                CJdp = population->dps[CJid-1];
                tempObj += (CIdp*CJdp) /
                                (sqrt(
                                      pow(abs(CIx + CIl / 2.0 - CJx - CJl / 2.0), 2)
                                    + pow(abs(CIy + CIw / 2.0 - CJy - CJw / 2.0), 2)
                                    + pow(abs(CIz + CIh / 2.0 - CJz - CJh / 2.0), 2)));
                if(couplings[CJid-1][CIid-1] == 1){
                    wireObj +=  abs(CIx + CIl / 2.0 - CJx - CJl / 2.0)
                            +   abs(CIy + CIw / 2.0 - CJy - CJw / 2.0)
                            +   (maxHeight + maxWidth) * abs(CIz + CIh / 2.0 - CJz - CJh / 2.0);
                }
            }
        }
    }
    if (isinf(tempObj)){
        tempObj=100000;
    }
    if (isinf(wireObj)){
        wireObj=100000;
    }
    res[0]= tempObj;
    res[1] = wireObj;
}

int Floorplanner::checkfeasibility(int indice_sol){
    int feasible=0;
    int totalSize = this->maxLength * this->maxWidth * this->maxHeight;
    int* occupied = new int[totalSize];
    for(int i=0;i<totalSize;i++){
        occupied[i] = 0;
    }
    for(int i=0;i<totalSize;i++){
        if(population->floorplan[indice_sol][i]>0){
            int idAux = population->floorplan[indice_sol][i];
            int lAux,wAux;
            if(!population->rotated[indice_sol][idAux-1]){
                lAux = population->ls[idAux-1];
                wAux = population->ws[idAux-1];
            }else{
                lAux = population->ws[idAux-1];
                wAux = population->ls[idAux-1];
            }
            
            //int hAux = solution->hs[idAux];
            int hAux = 1;
            int x,y,z;
            population->getCoordinates(&x,&y,&z,i);
            for(int l=0;l<lAux;l++){
                for(int w=0;w<wAux;w++){
                    for(int h=0;h<hAux;h++){
                        if(x+l>=maxLength)feasible++;
                        if(y+w>=maxWidth)feasible++;
                        if(z+h>=maxHeight)feasible++;
                        int index = i + w*maxLength + l + h*maxLength*maxWidth;
                        if((index<totalSize)&&(x+l<maxLength)&&(y+w<maxWidth)&&(z+h<maxHeight) ){
                            if(occupied[index]>0)feasible++;
                            occupied[index]++;
                        }else if(index>=totalSize){
                            feasible++;
                        }
                    }
                }
            }
        }
    }
    delete occupied;
    occupied = NULL;
    return feasible;
}

void Floorplanner::execute() {
    int nextPercentageReport = 10;
    while (currentGeneration < maxGenerations) {
        step();
        int percentage = round((currentGeneration * 100) / maxGenerations);
        if (percentage == nextPercentageReport) {
            cout << percentage << " % performed ...\t" ;//<< endl;
            cout << "feasibility: " << population->objectives[0][0] << endl;
            nextPercentageReport += 10;
        }
    }
}

void Floorplanner::step() {
    currentGeneration++;
    int ind_parent1;
    int ind_parent2;
    clearOffspring();
    for (int i = 0; i < (maxPopulation / 2); i++) {
        ind_parent1 = executeBinaryTournamentNSGAII();
        ind_parent2 = executeBinaryTournamentNSGAII();
        int ind_child1= 2*i + maxPopulation;
        int ind_child2= 2*i + 1 + maxPopulation;
        SinglePointCrossover(ind_parent1,ind_parent2,ind_child1,ind_child2);
        executeMutation(ind_child1);
        executeMutation(ind_child2);
        //evaluate(ind_child1);
        //evaluate(ind_child2);
    }
    evaluateChildren();
    population->rankSolutions();
    reduce(maxPopulation);

}

void Floorplanner::evaluateChildren(){
    int feasible[maxPopulation];
    for(int i=0;i<maxPopulation;i++){
        feasible[i] = checkfeasibility(i+maxPopulation);
    }
    
    float objTemp[maxPopulation];
    float objWire[maxPopulation];
    
    computeTempAndWireChildren(objTemp,objWire);

    for(int i=0;i<maxPopulation;i++){
        population->objectives[i+maxPopulation][0] = feasible[i];
        population->objectives[i+maxPopulation][1] = (feasible[i] + 1) * objTemp[i];
        population->objectives[i+maxPopulation][2] = (feasible[i] + 1) * objWire[i];
        if (feasible[i] == 0 && objTemp[i]<bestObjs[0] && objWire[i]<bestObjs[1]){
            cout << "One feasible solution found: (" << feasible[i] << "," << objTemp[i] << "," << objWire[i] << ")" << endl;
            bestObjs[0] = objTemp[i];
            bestObjs[1] = objWire[i];
        }
    }
}

void Floorplanner::clearOffspring(){
    int maxSize = maxPopulation;
    for(int i=maxSize;i<population->number_of_solutions;i++){
        for(int j=0;j<population->H*population->L*population->W;j++){
            population->floorplan[i][j]=0;
            population->surface[i][j]=0;
        }
        for(int j=0;j<population->number_of_elements;j++){
            //population->ls[i][j] = 0;
            //population->ws[i][j] = 0;
            population->rotated[i][j] = false;
            population->pos[i][j] = 0;
        }
        for(int j=0;j<population->number_of_objectives;j++){
            population->objectives[i][j]= 0;
        }
        population->crowdingDistances[i] = 0;
    }
}

int Floorplanner::executeBinaryTournamentNSGAII() {
    int popSize = population->number_of_solutions / 2;
    if (popSize < 2) {
        cout << "Population size must be greater or equal than 2." << endl;
        return -1;
    }
    int index1 = rand() % popSize;
    int index2 = index1;
    while (index2 == index1) {
        index2 = rand()%popSize;
    }
    int flag = population->compareSolutionDominance(index1,index2);
    if (flag < 0) {
        return index1;
    } else if (flag > 0) {
        return index2;
    } else if (population->crowdingDistances[index1] > population->crowdingDistances[index2]){
        return index1;
    } else if (population->crowdingDistances[index2] > population->crowdingDistances[index1]) {
        return index2;
    } else if ((rand() / (RAND_MAX + 1.0)) < 0.5 ) {
        return index1;
    } else {
        return index2;
    }
}

void Floorplanner::SinglePointCrossover(int ind_parent1,int ind_parent2,int ind_child1,int ind_child2) {
    int tam=population->number_of_elements;
    int indexPoint = rand() % tam;
    int maxL = population->L;
    int maxW = population->W;
    int maxH = population->H;
    int totalSize=maxL*maxW*maxH;
    for(int i=0; i<indexPoint;++i){
        int id = population->ids[i];
        int index_corners = population->GetIndexFromId(ind_parent1,id);
        while(population->floorplan[ind_child1][index_corners] > 0 || (!population->feasiblePosition(ind_child1,id,index_corners))){
            index_corners = rand() % (totalSize);
        }
        population->floorplan[ind_child1][index_corners] = id;
        population->pos[ind_child1][id-1] = index_corners;
        //population->ls[ind_child1][id-1] = population->ls[ind_parent1][id-1];
        //population->ws[ind_child1][id-1] = population->ws[ind_parent1][id-1];
        population->rotated[ind_child1][id-1] = population->rotated[ind_parent1][id-1];
        index_corners = population->GetIndexFromId(ind_parent2,id);
        while(population->floorplan[ind_child2][index_corners] > 0  || (!population->feasiblePosition(ind_child2,id,index_corners))){
            index_corners = rand() % (totalSize);
        }
        population->floorplan[ind_child2][index_corners] = id;
        population->pos[ind_child2][id-1] = index_corners;
        //population->ls[ind_child2][id-1] = population->ls[ind_parent2][id-1];
        //population->ws[ind_child2][id-1] = population->ws[ind_parent2][id-1];
        population->rotated[ind_child2][id-1] = population->rotated[ind_parent2][id-1];
    }

    for(int i=indexPoint;i<tam;i++){
        int id = population->ids[i];
        int index_corners = population->GetIndexFromId(ind_parent2,id);
        while(population->floorplan[ind_child1][index_corners] > 0 || (!population->feasiblePosition(ind_child1,id,index_corners))){
            index_corners = rand() % (totalSize);
        }
        population->floorplan[ind_child1][index_corners] = id;
        population->pos[ind_child1][id-1] = index_corners;
        //population->ls[ind_child1][id-1] = population->ls[ind_parent2][id-1];
        //population->ws[ind_child1][id-1] = population->ws[ind_parent2][id-1];
        population->rotated[ind_child1][id-1] = population->rotated[ind_parent2][id-1];
        index_corners = population->GetIndexFromId(ind_parent1,id);
        while(population->floorplan[ind_child2][index_corners] > 0  || (!population->feasiblePosition(ind_child2,id,index_corners))){
            index_corners = rand() % (totalSize);
        }
        population->floorplan[ind_child2][index_corners] = id;
        population->pos[ind_child2][id-1] = index_corners;
        //population->ls[ind_child2][id-1] = population->ls[ind_parent1][id-1];
        //population->ws[ind_child2][id-1] = population->ws[ind_parent1][id-1];
        population->rotated[ind_child2][id-1] = population->rotated[ind_parent1][id-1];
    }
}

void Floorplanner::executeMutation(int solInd) {
    int oldPos,newPos;
    float random;
    int totalSize = population->H * population->W * population->L;
    for(int i=0;i<population->number_of_elements;i++){
        int Iid = population->ids[i];
        oldPos=population->pos[solInd][Iid-1];
        newPos=0;
        if ((rand() / (RAND_MAX + 1.0)) < this->mutation_probability){
            random = rand()/(RAND_MAX + 1.0);
            if(random < 0.5){
                newPos = rand() % totalSize;
                while(!population->feasiblePosition(solInd,Iid,newPos)){
                    newPos = rand() % totalSize;
                }
                if(population->floorplan[solInd][newPos] > 0){
                    int Jid=population->floorplan[solInd][newPos];
                    population->pos[solInd][Iid-1]=newPos;
                    population->floorplan[solInd][newPos] = Iid;
                    while(population->floorplan[solInd][oldPos]>0 || (!population->feasiblePosition(solInd,Jid,oldPos))){
                        oldPos = rand() % totalSize;
                    }
                    population->floorplan[solInd][oldPos] = Jid;
                    population->pos[solInd][Jid-1] = oldPos;
                } else{
                    population->floorplan[solInd][oldPos] = 0;
                    population->floorplan[solInd][newPos] = Iid;
                    population->pos[solInd][Iid-1] = newPos;
                }
            } else if (random < 1.00) {
                int l,w;
                if(!population->rotated[solInd][Iid-1]){
                    l = population->ls[Iid-1];
                    w = population->ws[Iid-1];
                }else{
                    l = population->ws[Iid-1];
                    w = population->ls[Iid-1];
                }
                int h=1;
                int x,y,z;
                int index = population->GetIndexFromId(solInd,Iid);
                population->getCoordinates(&x,&y,&z,index);
                int feasible = 0;
                if(x+w > population->L)feasible++;
                if(y+l > population->W)feasible++;
                if(z+h > population->H)feasible++;
                if(feasible==0){
                    //population->ls[solInd][Iid-1] = w;
                    //population->ws[solInd][Iid-1] = l;
                    bool bAux = population->rotated[solInd][Iid-1];
                    population->rotated[solInd][Iid-1] = (!bAux);
                }
            }
        }
    }
}

void Floorplanner::reduce(int maxSize) {
    int indexInsert=0;
    int indRank=0;
    while(indexInsert<maxSize){
        int lookforRank = population->rankValues[indRank];
        int j=indexInsert;
        for(int i=j;i<population->number_of_solutions;i++){
            if(population->ranks[i]==lookforRank){
                population->replace(indexInsert,i,indRank);
                indexInsert++;
            }
        }
        indRank++;
    }
    for(int i=0;i<population->number_of_solutions;i++){
        population->ranks[i] = population->ranks[i] + 1;
    }
    for(int i=maxSize;i<population->number_of_solutions;i++){
        for(int j=0;j<population->H*population->L*population->W;j++){
            population->floorplan[i][j]=0;
        }
        for(int j=0;j<population->number_of_elements;j++){
            //population->ls[i][j] = 0;
            //population->ws[i][j] = 0;
            population->rotated[i][j] = false;
            population->pos[i][j] = 0;
        }
        for(int j=0;j<population->number_of_objectives;j++){
            population->objectives[i][j]= 0;
        }
        population->crowdingDistances[i] = 0;
    }
}

void Floorplanner::saveNonDominatedFront(){
    for(int i=0;i<maxPopulation;i++){
        if(population->ranks[i]==1){
            printSolution(i);
            population->solutionToXml(i);
        }
    }
}

void Floorplanner::printSolution(int indSol){
    /*int L = population->L;
    int W = population->W;
    int H = population->H;*/
    cout << "SOLUTION " << indSol << endl;
    //cout << "Dimensions: L=" << L << ", W=" << W << ", H=" << H << endl;
    cout << "Objectives: FEASIBILITY: " <<population->objectives[indSol][0] << ",TEMP: " << population->objectives[indSol][1] << ", WIRE: "<< population->objectives[indSol][2] << " " << endl;
    checkfeasibility(indSol);
    printFloorplan(indSol);
    cout << endl;
}

void Floorplanner::printFloorplan(int indSol){
    int feasible=0;
    int L = population->L;
    int W = population->W;
    int H = population->H;
    int totalSize = L * W * H;
    int* occupied = new int[totalSize];
    for(int i=0;i<totalSize;i++){
        occupied[i] = 0;
    }
    for(int i=0;i<totalSize;i++){
        if(population->floorplan[indSol][i]>0){
            int idAux = population->floorplan[indSol][i];
            int lAux,wAux;
            if(!population->rotated[indSol][idAux-1]){
                lAux = population->ls[idAux-1];
                wAux = population->ws[idAux-1];
            }else{
                lAux = population->ws[idAux-1];
                wAux = population->ls[idAux-1];
            }
            //int hAux = solution->hs[idAux];
            int hAux = 1;
            int x,y,z;
            population->getCoordinates(&x,&y,&z,i);
            for(int l=0;l<lAux;l++){
                for(int w=0;w<wAux;w++){
                    for(int h=0;h<hAux;h++){
                        if(x+l>=L)feasible++;
                        if(y+w>=W)feasible++;
                        //if(z+h>H)feasible++;
                        int index = i + w*L + l + h*L*W;
                        if((index<totalSize)&&(x+l<L)&&(y+w<W)&&(z+h<H) ){
                            if(occupied[index]>0){
                                occupied[index] = -1;
                            }else if(occupied[index]!=-1){
                                occupied[index] = idAux;
                            }
                        }else if(index>=totalSize){
                            feasible++;
                        }
                    }
                }
            }
        }
    }

    for(int i=1; i<L*W*H+1;i++){
        printf("%02d ", occupied[i-1]);
        if(i % L == 0) cout << endl;
        if(i % (W*L) ==0) cout << endl;
    }
    cout << "Out of the borders: " << feasible << endl;
    delete occupied;
    occupied = NULL;
}

//int main(int argc, char** argv) {
int main() {
    //srand(time(NULL));
    srand(2001);
    cout << "Parameters: fileInput numLayers NumInd NumGen" << endl; 
    string fileInput = "test/dp_sumaPond_30.xml";
    int maxPopulation = 128;
    int maxGenerations = 1;
    int numObjectives = 3;
    Floorplanner* problem = new Floorplanner(fileInput,maxPopulation,maxGenerations, numObjectives);
    ComponentVariable* cv=problem->components.at(0).clone();
    problem->initialize();
    problem->execute();
    //problem->saveNonDominatedFront();
    delete cv; cv=NULL;
    delete problem; problem=NULL;
    exit(0);
}
