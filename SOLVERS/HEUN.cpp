#include "../SOLVERS/HEUN.h"
#include "../SIMTIME.h"
#include "../../CMAT/Vector.h"
#include "../SSMODEL.h"

Heun::Heun(StateSpaceModel *model) : Solver(model)
{
}

Heun::~Heun()
{

}

void Heun::integrate(SimTime &st)
{
    st.setReal(false);

    double t = st.t;
    double h = st.dt;

    Vector x = model->getContinuousStatesValueVector();
    Vector k1 = h * model->getContinuousStatesDerivativeVector();

    st.t = (st.cnt + 1) * h;
    model->setContinuousStatesValueVector(x + k1);
    model->output(st);
    Vector k2 = h * model->getContinuousStatesDerivativeVector();

    model->setContinuousStatesValueVector(x + (k1 + k2) / 2);

    st.t = t;
}
