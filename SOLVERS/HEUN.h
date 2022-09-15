#ifndef HEUN_H
#define HEUN_H

#include "../SOLVER.h"

class Heun : public Solver
{
public:
    Heun(StateSpaceModel* model);
    virtual ~Heun();

protected:
    virtual void integrate(SimTime& st) = 0;
};

#endif // HEUN_H
