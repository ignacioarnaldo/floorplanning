#include "Floorplanner.h"
#include <shrUtils.h>
#include "cutil_inline.h"

#include <computeWireTemp_kernel.cu>

/*void Floorplanner::computeTempAndWireChildren(float* objTemp, float* objWire){
    float tempObj = 0;
    float wireObj = 0;
    int CIz,CIy,CIx;
    int CJz,CJy,CJx;
    int CIid,posI,CIl,CIw,CIh;
    float CIdp;
    int CJid,posJ,CJl,CJw,CJh;
    float CJdp;
    for(int s=maxPopulation;s<2*maxPopulation;s++){
	tempObj = 0;
	wireObj = 0;
        for (int i = 0; i < numberOfVariables; ++i) {
            CIid = population->ids[i];
            //posI = population->GetIndexFromId(s,CIid);
	    posI = population->pos[s][CIid-1];
            population->getCoordinates(&CIx,&CIy,&CIz,posI);
            CIl = population->ls[s][CIid-1];
            CIw = population->ws[s][CIid-1];
            //int CIh = solutions->hs[indice_sol][i];
            CIh = 1;
            CIdp = population->dps[CIid-1];
            for (int j = 0; j < numberOfVariables; ++j) {
                if (i < j) {
                    CJid = population->ids[j];
                    //posJ = population->GetIndexFromId(s,CJid);
		    posJ = population->pos[s][CJid-1];
                    population->getCoordinates(&CJx,&CJy,&CJz,posJ);
                    CJl = population->ls[s][CJid-1];
                    CJw = population->ws[s][CJid-1];
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
        objTemp[s-maxPopulation] = tempObj;
        objWire[s-maxPopulation] = wireObj;
    }
}*/

void Floorplanner::computeTempAndWireChildren(float* objTemp, float* objWire){

    int mem_size_objTemp = sizeof(float) * maxPopulation;
    int mem_size_objWire = sizeof(float) * maxPopulation;
    int mem_size_ids = sizeof(int) * numberOfVariables;
    int mem_size_dps = sizeof(float) * numberOfVariables;
    int mem_size_ls = sizeof(int) * numberOfVariables;
    int mem_size_ws = sizeof(int) * numberOfVariables;
    int mem_size_rotated = sizeof(bool) * maxPopulation * numberOfVariables;
    int mem_size_pos = sizeof(int) * maxPopulation * numberOfVariables;
    int mem_size_couplings = sizeof(bool) * numberOfVariables * numberOfVariables;

    int* h_ls = new int[numberOfVariables];
    int* h_ws = new int[numberOfVariables];
    bool* h_rotated = new bool[maxPopulation * numberOfVariables];
    int* h_pos = new int[maxPopulation * numberOfVariables];
    int* h_couplings = new int[numberOfVariables * numberOfVariables];
    
    for(int j=0;j<numberOfVariables;j++){
        h_ls[j] = population->ls[j];
	h_ws[j] = population->ws[j];
    }

    for(int i=0;i<maxPopulation;i++){
	for(int j=0;j<numberOfVariables;j++){
	    //h_ls[i*numberOfVariables + j] = population->ls[i+maxPopulation][j];
	    //h_ws[i*numberOfVariables + j] = population->ws[i+maxPopulation][j];
	    h_rotated[i*numberOfVariables + j] = population->rotated[i+maxPopulation][j];
	    h_pos[i*numberOfVariables + j] = population->pos[i+maxPopulation][j];
	}
    }
    for(int i=0;i<numberOfVariables;i++){
	for(int j=0;j<numberOfVariables;j++){
	    h_couplings[i*numberOfVariables + j] = couplings[i][j];
	}
    }
    
    float* d_objTemp;
    float* d_objWire;
    int* d_ids;
    float* d_dps;
    int* d_ls;
    int* d_ws;
    bool* d_rotated;
    int* d_pos;
    int* d_couplings;
	
    // allocate device memory
    cutilSafeCall(cudaMalloc((void**) &d_objTemp, mem_size_objTemp));
    cutilSafeCall(cudaMalloc((void**) &d_objWire, mem_size_objWire));
    cutilSafeCall(cudaMalloc((void**) &d_ids, mem_size_ids));
    cutilSafeCall(cudaMalloc((void**) &d_dps, mem_size_dps));
    cutilSafeCall(cudaMalloc((void**) &d_ls, mem_size_ls));
    cutilSafeCall(cudaMalloc((void**) &d_ws, mem_size_ws));
    cutilSafeCall(cudaMalloc((void**) &d_rotated, mem_size_rotated));
    cutilSafeCall(cudaMalloc((void**) &d_pos, mem_size_pos));
    cutilSafeCall(cudaMalloc((void**) &d_couplings, mem_size_couplings));

    //set Timer Transfer input
    unsigned int timerInput = 0;
    cutilCheckError(cutCreateTimer(&timerInput));
    cutilCheckError(cutStartTimer(timerInput));

    // copy host memory to device
    cutilSafeCall(cudaMemcpy(d_ids, population->ids, mem_size_ids,cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpy(d_dps, population->dps, mem_size_dps,cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpy(d_ls,  h_ls,  mem_size_ls, cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpy(d_ws,  h_ws,  mem_size_ws, cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpy(d_rotated, h_rotated,  mem_size_rotated, cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpy(d_pos, h_pos, mem_size_pos, cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpy(d_couplings, h_couplings, mem_size_couplings, cudaMemcpyHostToDevice) );

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timerInput));
    double dSecondsInput = cutGetTimerValue(timerInput);
    //cout << "Transfer input time: " << dSecondsInput << endl;

    // setup execution parameters
    int BLOCK_SIZE_X = 128;
    int BLOCK_SIZE_Y = 1;
    dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    int gridx = 1;
    int gridy = 1;
    dim3 grid(gridx, gridy);

    //set Timer Execution
    unsigned int timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));


    //launch kernel    
    computeWireTempKernel<<< grid, threads >>>(maxPopulation,numberOfVariables,d_objTemp,d_objWire,d_ids,
					       d_dps,d_ls,d_ws,d_rotated,d_pos,d_couplings,maxHeight,maxLength,maxWidth);

    /*computeWireTempKernelCPU(maxPopulation,numberOfVariables,objTemp,objWire,population->ids,population->dps,population->ls,
			      population->ws,population->pos,this->couplings,maxHeight,maxLength,maxWidth);*/

    //thread synchronization
    cudaThreadSynchronize();

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    double dSeconds = cutGetTimerValue(timer);
    //cout << "Execution time: " << dSeconds << endl;
 
    //set Timer Transfer output
    unsigned int timerOutput = 0;
    cutilCheckError(cutCreateTimer(&timerOutput));
    cutilCheckError(cutStartTimer(timerOutput));
 
    // copy device to host memory
    cutilSafeCall(cudaMemcpy(objTemp, d_objTemp, mem_size_objTemp,cudaMemcpyDeviceToHost) );
    cutilSafeCall(cudaMemcpy(objWire, d_objWire, mem_size_objWire,cudaMemcpyDeviceToHost) );

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timerOutput));
    double dSecondsOutput = cutGetTimerValue(timerOutput);
    //cout << "Transfer output time: " << dSecondsOutput << endl;
	
    cout << "Max Population: " << maxPopulation << endl;
    for(int i=0;i<maxPopulation;i++){
	cout << "Individual " << i << " objs: " << objTemp[i] << ", " << objWire[i] << endl;
    }

    //free device memory
    cutilSafeCall(cudaFree(d_objTemp));
    cutilSafeCall(cudaFree(d_objWire));
    cutilSafeCall(cudaFree(d_ids));
    cutilSafeCall(cudaFree(d_dps));
    cutilSafeCall(cudaFree(d_ls));
    cutilSafeCall(cudaFree(d_ws));
    cutilSafeCall(cudaFree(d_pos));
    cutilSafeCall(cudaFree(d_couplings));

    //free host memory
    delete[] h_ls;
    delete[] h_ws;
    delete[] h_pos;
    delete[] h_couplings;
    cudaThreadExit();
}

void Floorplanner::computeWireTempKernelCPU(int maxPop,int numberOfVars,float* objT, float* objW,int* ids,float* dps,int* ls,
			      int* ws,bool** rotated,int** pos,bool** couplings,int H,int L,int W){
    int CIz,CIy,CIx;
    int CJz,CJy,CJx;
    int CIid,posI,CIl,CIw,CIh;
    float CIdp;
    int CJid,posJ,CJl,CJw,CJh;
    float CJdp;
    float tempObj = 0;
    float wireObj = 0;
    for(int s=maxPop;s<2*maxPop;s++){
	tempObj = 0;
	wireObj = 0;
        for (int i = 0; i < numberOfVars; ++i) {
            CIid = ids[i];
	    posI = pos[s][CIid-1];
            //getCoordinates(&CIx,&CIy,&CIz,posI);
	    int totalsize = L * W * H;
	    int cellsPerLayer = totalsize / H;
	    CIz = (posI / cellsPerLayer);
	    posI = posI - (cellsPerLayer * CIz);
	    CIx = posI % L;
	    posI = posI - CIx;
	    CIy = posI / L;
            if(!rotated[s][CIid-1]){
                CIl = ls[CIid-1];
                CIw = ws[CIid-1];
            }else{
                CIl = ws[CIid-1];
                CIw = ls[CIid-1];
            }
            //int CIh = solutions->hs[indice_sol][i];
            CIh = 1;
            CIdp = dps[CIid-1];
            for (int j = 0; j < numberOfVars; ++j) {
                if (i < j) {
                    CJid = ids[j];
  		    posJ = pos[s][CJid-1];
                    //getCoordinates(&CJx,&CJy,&CJz,posJ);
		    //int totalsize = L * W * H;
		    //int cellsPerLayer = totalsize / H;
		    CJz = (posJ / cellsPerLayer);
		    posJ = posJ - (cellsPerLayer * CJz);
		    CJx = posJ % L;
		    posJ = posJ - CJx;
		    CJy = posJ / L;
                    if(!rotated[s][CJid-1]){
                        CJl = ls[CJid-1];
                        CJw = ws[CJid-1];
                    }else{
                        CJl = ws[CJid-1];
                        CJw = ls[CJid-1];
                    }
                    //int CJh = solution->hs[j];
                    CJh = 1;
                    CJdp = dps[CJid-1];
                    tempObj += (CIdp*CJdp) /
                                    (sqrt(
                                          pow(abs(CIx + CIl / 2.0 - CJx - CJl / 2.0), 2)
                                        + pow(abs(CIy + CIw / 2.0 - CJy - CJw / 2.0), 2)
                                        + pow(abs(CIz + CIh / 2.0 - CJz - CJh / 2.0), 2)));
                    if(couplings[CJid-1][CIid-1] == 1){
                        wireObj +=  abs(CIx + CIl / 2.0 - CJx - CJl / 2.0)
                                +   abs(CIy + CIw / 2.0 - CJy - CJw / 2.0)
                                +   (H + W) * abs(CIz + CIh / 2.0 - CJz - CJh / 2.0);
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
        objT[s-maxPop] = tempObj;
        objW[s-maxPop] = wireObj;
    }
}
