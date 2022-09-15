#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

#include "../SOLVER.h"

class RungeKutta4 : public Solver
{
public:
    RungeKutta4(StateSpaceModel* model);
    virtual ~RungeKutta4();

protected:
    void integrate(SimTime& st);
};

#endif // RUNGEKUTTA4_H
