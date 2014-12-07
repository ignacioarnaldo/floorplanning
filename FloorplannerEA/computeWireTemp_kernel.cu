#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
//! Matrix multiplication on the device: LHS = A * V
//! wA is A's width and wV is V's width
////////////////////////////////////////////////////////////////////////////////

__global__ void computeWireTempKernel(int maxPop,int numberOfVars,float* objT, float* objW,int* ids,float* dps,int* ls,
				      int* ws,bool* rotated,int* pos,int* couplings,int H,int L,int W){

	int i=blockIdx.x*blockDim.x+threadIdx.x;
	int j=blockIdx.y*blockDim.y+threadIdx.y;
	int index = i + j * blockDim.x;
	__shared__ int sids[61];
	__shared__ int sls[61];
	__shared__ int sws[61];
	__shared__ float sdps[61];
	__shared__ int spos[128*61];
	__shared__ bool scouplings[128*61];
	__shared__ bool srotated[128*61];
	if(index<61){
		sids[index] = ids[index];
		//__syncthreads();
	}
__syncthreads();
	if(index<61){
		sls[index] = ls[index];
		//__syncthreads();
	}
__syncthreads();
	if(index<61){
		sws[index] = ws[index];
		//__syncthreads();
	}
__syncthreads();
	if(index<61){
		sdps[index] = dps[index];
		//__syncthreads();
	}
__syncthreads();	
	for(int i=0;i<61;i++){
		spos[i*61 + index] = pos[i*61 + index];
		//__syncthreads();
	}
__syncthreads();
	for(int i=0;i<61;i++){
		scouplings[i*61 + index] = couplings[i*61 + index];
		//__syncthreads();
	}
__syncthreads();
	for(int i=0;i<61;i++){
		srotated[i*61 + index] = rotated[i*61 + index];
		//__syncthreads();
	}
__syncthreads();
	float tempObj = 0;
	float wireObj = 0;
	int CIz,CIy,CIx,CJz,CJy,CJx,CIid,posI,CIl,CIw,CIh,CJid,posJ,CJl,CJw,CJh;
        float CIdp,CJdp;
        int totalsize = L * W * H;
	int cellsPerLayer = totalsize / H;
        for (int i = 0; i < numberOfVars; ++i) {
            CIid = sids[i];
            CIdp = sdps[CIid-1];
	    posI = spos[index*numberOfVars + CIid-1];
	    if(!srotated[index*numberOfVars + CIid-1]){
	        CIl = sls[CIid-1];
        	CIw = sws[CIid-1];
	    }else{
		CIl = sws[CIid-1];
        	CIw = sls[CIid-1];
	    }
	    CIz = (posI / cellsPerLayer);
	    posI = posI - (cellsPerLayer * CIz);
	    CIx = posI % L;
	    posI = posI - CIx;
	    CIy = posI / L;
            CIh = 1;
            for (int j = 0; j < numberOfVars; ++j) {
                if (i < j) {
                    CJid = sids[j];
		    CJdp = sdps[CJid-1];		    
		    posJ = spos[index*numberOfVars + CJid-1];
		    if(!srotated[index*numberOfVars + CJid-1]){
		        CJl = sls[CJid-1];
	        	CJw = sws[CJid-1];
		    }else{
			CJl = sls[CJid-1];
	        	CJw = sws[CJid-1];
		    }		    
		    CJz = (posJ / cellsPerLayer);
		    posJ = posJ - (cellsPerLayer * CJz);
		    CJx = posJ % L;
		    posJ = posJ - CJx;
		    CJy = posJ / L;
                    CJh = 1;
                    
                    tempObj += (CIdp*CJdp) /
                                    (sqrtf(
                                          pow(abs(CIx + CIl / 2.0 - CJx - CJl / 2.0), 2)
                                        + pow(abs(CIy + CIw / 2.0 - CJy - CJw / 2.0), 2)
                                        + pow(abs(CIz + CIh / 2.0 - CJz - CJh / 2.0), 2)));
		    
                    if(scouplings[(CJid-1)*numberOfVars + CIid-1] == true){
		       wireObj += (fabsf(CIx+CIl/2.0-CJx-CJl/2.0) + fabsf(CIy+CIw/2.0-CJy-CJw/2.0) + (H+W)*fabsf(CIz+CIh/2.0-CJz-CJh/2.0));
  		    }
                }
		
            }
        }
	//__syncthreads();
	objT[index] = tempObj;
	objW[index] = index;	
}

