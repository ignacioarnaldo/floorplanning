#include "ComponentVariable.h"

ComponentVariable::ComponentVariable(int id, string name, int type, int x, int y, int z, int l, int w, int h, float dp, int cellSize) {
        this->id = id;
        this->name = name;
        this->type = type;
        this->x = x;
        this->y = y;
        this->z = z;
        this->l = l;
        this->w = w;
        this->h = h;
        //this->p = p;
        //this->dp = p/(l*w*h);
        this->dp = dp;
        this->cellSize = cellSize;
        this->p = dp * l * w * 0.000001 * 0.000001 * cellSize * cellSize;
}

ComponentVariable* ComponentVariable::clone() {
    ComponentVariable* clone = new ComponentVariable(this->id, this->name, this->type, this->x, this->y, this->z, this->l, this->w, this->h, this->dp,this->cellSize);
    return clone;
}
