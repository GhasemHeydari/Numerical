#include "RK4.h"
#include "../SIMTIME.h"
#include "../../CMAT/Vector.h"
#include "../SSMODEL.h"
#include <stdio.h>

RungeKutta4::RungeKutta4(StateSpaceModel *model) : Solver(model)
{
}

RungeKutta4::~RungeKutta4()
{

}

void RungeKutta4::integrate(SimTime &st) {
    st.setReal(false);

    double t = st.t;
    double h = st.dt;

    Vector x = model->getContinuousStatesValueVector();
    Vector k1 = h * model->getContinuousStatesDerivativeVector();
    st.t = t + h / 2;
    model->setContinuousStatesValueVector(x + k1 / 2);
    model->output(st);
    Vector k2 = h * model->getContinuousStatesDerivativeVector();

    st.t = t + h / 2;
    model->setContinuousStatesValueVector(x + k2 / 2);
    model->output(st);
    Vector k3 = h * model->getContinuousStatesDerivativeVector();

    st.t = (st.cnt + 1) * h;
    model->setContinuousStatesValueVector(x + k3);
    model->output(st);
    Vector k4 = h * model->getContinuousStatesDerivativeVector();

    model->setContinuousStatesValueVector(x + (k1 + 2 * k2 + 2 * k3 + k4 ) / 6);

    st.t = t;
}
