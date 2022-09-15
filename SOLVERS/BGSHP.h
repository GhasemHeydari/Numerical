#ifndef BOGACKISHAMPINE_H
#define BOGACKISHAMPINE_H

#include "../SOLVER.h"

class BogackiShampine : public Solver
{
public:
    BogackiShampine(StateSpaceModel* model);
    virtual ~BogackiShampine();

protected:
    void integrate(SimTime& st);
};

#endif // BOGACKISHAMPINE_H
