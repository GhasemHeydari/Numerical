#ifndef EULER_H
#define EULER_H

#include "../SOLVER.h"
#include "../../CMAT/Vector.h"
#include "../SSMODEL.h"
#include "../SIMTIME.h"
#include "../CSTATE.h"

class Euler : public Solver
{
public:
    Euler(StateSpaceModel* model);
    virtual ~Euler();

protected:
    void integrate(SimTime& st);
};

#endif // EULER_H
