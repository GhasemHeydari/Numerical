#include "EULER.h"

Euler::Euler(StateSpaceModel *model) : Solver(model)
{
}

Euler::~Euler()
{

}

void Euler::integrate(SimTime &st) {
    for(int i = 0; i < model->continuousScalarStatesCount; i++) {
        ContinuousState* cs = model->continuousScalarStatesList[i];
        cs->setValue(cs->getValue() + st.dt * cs->getDerivative());
    }
}
