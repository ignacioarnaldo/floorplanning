#include <sstream>

using namespace std;

#ifndef COMPONENTVARIABLE_H
#define	COMPONENTVARIABLE_H

class ComponentVariable {
public:
    ComponentVariable() {}
    ~ComponentVariable() {name="";}
    ComponentVariable(int id, string name, int type, int x, int y, int z, int l, int w, int h, float p,int cellSize);
    ComponentVariable* clone();

    int id;
    string name;
    int type; // 0->Core, 1->Rest
    int x;
    int y;
    int z;
    int l;
    int w;
    int h;
    float p;
    float dp; // Power density
    int cellSize;
protected:
    
};

#endif	/* COMPONENTVARIABLE_H */


