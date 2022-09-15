#include "BGSHP.h"
#include "../../CMAT/Vector.h"
#include "../SIMTIME.h"
#include "../SSMODEL.h"

const double DIV_1_2 = 1.0 / 2;
const double DIV_3_4 = 3.0 / 4;
const double DIV_2_9 = 2.0 / 9;
const double DIV_1_3 = 1.0 / 3;
const double DIV_4_9 = 4.0 / 9;

BogackiShampine::BogackiShampine(StateSpaceModel *model) : Solver(model)
{
}

BogackiShampine::~BogackiShampine()
{

}

void BogackiShampine::integrate(SimTime &st) {
    st.setReal(false);

    double t = st.t;
    double h = st.dt;

    Vector x = model->getContinuousStatesValueVector();
    Vector k1 = h * model->getContinuousStatesDerivativeVector();

    st.t = t + DIV_1_2 * h;
    model->setContinuousStatesValueVector(x + DIV_1_2 * k1);
    model->output(st);
    Vector k2 = h * model->getContinuousStatesDerivativeVector();

    st.t = t + DIV_3_4 * h;
    model->setContinuousStatesValueVector(x + DIV_3_4 * k2);
    model->output(st);
    Vector k3 = h * model->getContinuousStatesDerivativeVector();

    model->setContinuousStatesValueVector(x + DIV_2_9 * k1 + DIV_1_3 * k2 + DIV_4_9 * k3);

    st.t = t;
}
