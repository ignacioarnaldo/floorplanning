#include "FloorplanConfiguration.h"


FloorplanConfiguration::FloorplanConfiguration() {
}

FloorplanConfiguration::~FloorplanConfiguration() {
    
}

FloorplanConfiguration::FloorplanConfiguration(string filePath) {
    this->filePath = filePath;
    components.clear();
    CMarkup cmark;
    cmark.Load(filePath);
    cmark.ResetPos();
    cmark.FindElem();
    cellSize = string2int(cmark.GetAttrib("CellSize"));
    maxLength = string2int(cmark.GetAttrib("Length"));
    maxWidth = string2int(cmark.GetAttrib("Width"));
    maxHeight = string2int(cmark.GetAttrib("NumLayers"));
    numPowerProfiles = string2int(cmark.GetAttrib("NumPowerProfiles"));
    cmark.IntoElem();
    cmark.FindElem("Blocks");
    cmark.IntoElem();
    int area_occupied=0;
    int max_id=0;
    while (cmark.FindElem("Block")){
        int id = string2int(cmark.GetAttrib("id"));
        string name = cmark.GetAttrib("name");
        int type = string2int(cmark.GetAttrib("type"));
        string turnAsString = cmark.GetAttrib("turn");
        int turn=-1;
        if(turnAsString.length()>0) turn = string2int(turnAsString);
        int x = string2int(cmark.GetAttrib("x"));
        int y = string2int(cmark.GetAttrib("y"));
        int z = string2int(cmark.GetAttrib("z"));
        int l = string2int(cmark.GetAttrib("l"));
        int w = string2int(cmark.GetAttrib("w"));
        int h = string2int(cmark.GetAttrib("h"));
        //float p = string2float(cmark.GetAttrib("p"));
        float dp[numPowerProfiles];
        stringstream ss;
        for(int i=0;i<numPowerProfiles;i++){
            ss << "dp" << i;
            //dp[i] = string2float(cmark.GetAttrib("dp"+i));
            dp[i] = string2float(cmark.GetAttrib(ss.str()));
        }
        ComponentVariable component = ComponentVariable(id, name, type, x, y, z, l, w, h, dp[0],cellSize);
        components.insert(pair<int,ComponentVariable>(component.id,component));
        area_occupied += l * w * h;
        if(id>max_id)max_id=id;
    }

    numComponents = max_id;
    couplings = new bool* [numComponents];
    couplings[0] = new bool[numComponents*numComponents];
    for (int i = 1; i < numComponents; ++i){
        couplings[i] = couplings[i-1] + numComponents;
    }
    for(int i=0;i<numComponents;i++){
        for(int j=0;j<numComponents;j++){
            couplings[i][j]=false;
        }
    }

    float maxDP = DBL_MIN;
    set<int> keys;
    getKeySetMapComp(&components,&keys);
    set<int>::const_iterator iter;
    for (iter = keys.begin(); iter != keys.end(); ++iter){
        ComponentVariable component= components.at(*iter);
        if(maxDP<component.dp) maxDP = component.dp;
    }
    for (iter = keys.begin(); iter != keys.end(); ++iter){
        ComponentVariable component= components.at(*iter);
        float dpAux=component.dp;
        component.dp = (dpAux/maxDP);
        components.at(*iter)=component;
    }
    cmark.OutOfElem();
    cmark.FindElem("Couplings");
    cmark.IntoElem();
    while (cmark.FindElem("Coupling")){
        int idFrom = string2int(cmark.GetAttrib("idFrom"));
        int idTo = string2int(cmark.GetAttrib("idTo"));
        couplings[idFrom-1][idTo-1]=1;
        couplings[idTo-1][idFrom-1]=1;
    }
}

void FloorplanConfiguration::getKeySetMapComp(map<int,ComponentVariable>* mapcomps,set<int>* keys){
    keys->clear();
    map<int,ComponentVariable>::const_iterator itr;
    itr=mapcomps->begin();
    while (mapcomps->end() != itr){
        keys->insert((itr++)->first);
    }
}

int FloorplanConfiguration::string2int(string strConvert) {
  int intReturn;
  intReturn = atoi(strConvert.c_str());
  return(intReturn);
}

float FloorplanConfiguration::string2float(string strConvert){
    float floatReturn;
    floatReturn = atof(strConvert.c_str());
    return floatReturn;
}