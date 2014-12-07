#include <map>
#include <set>
#include "ComponentVariable.h"
#include <cfloat>
#include "Markup.h"

using namespace std;

#ifndef FLOORPLANCONFIGURATION_H
#define	FLOORPLANCONFIGURATION_H

class FloorplanConfiguration {
public:
    FloorplanConfiguration();
    FloorplanConfiguration(string filePath);
    ~FloorplanConfiguration();
    int string2int(string strConvert);
    float string2float(string strConvert);
    void getKeySetMapComp(map<int,ComponentVariable>* mapcomps,set<int>* keys);
    
    string filePath;
    int cellSize, maxLength, maxWidth, maxHeight,numPowerProfiles,numComponents;
    map<int,ComponentVariable> components;
    bool** couplings;
};

#endif	/* FLOORPLANCONFIGURATION_H */