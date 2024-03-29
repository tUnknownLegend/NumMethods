#include <cmath>
#include <iostream>
#include "shared.h"
#include "nonLinearSolve.h"
#include "SOLE.h"

using namespace std;


enum calcMethod {
    MexplicitEuler, MimplicitEuler, Msymmetric, MrungeKutta2, MrungeKutta4, MrungeKuttaRungeStep
};

vector<vector<TT>> calcDiff(const vector<vector<TT>> &answer) {
    vector<vector<TT>> diff(answer);

    for (int i = 0; i < answer.size(); ++i) {
        diff[i][0] = abs(cos(i * step) - answer[i][0]);
        diff[i][1] = abs(-sin(i * step) - answer[i][1]);
    }
    return diff;
}

vector<vector<TT>> explicitEuler(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; ++i) {
        y[i + 1] = vectorOperation(y[i],
                                   vectorRDigit(step, f(y[i]), '*'), '+');
    }
    return y;
}

vector<vector<TT>> implicitEuler(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = Newton("Euler", cond.size(), y[i]);
    }
    return y;
}

vector<vector<TT>> symmetric(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = Newton("Symmetric", cond.size(), y[i]);
    }
    return y;
}

void k_iSum(vector<TT> &temp, const vector<TT> &k1, const vector<TT> &k2, const vector<TT> &k3,
            const vector<TT> &k4, const vector<TT> &y_i, const TT Tchange) {
    temp = vectorOperation(vectorRDigit(2.0, k2, '*'), vectorRDigit(2.0, k3, '*'), '+');
    temp = vectorOperation(k1, temp, '+');
    temp = vectorOperation(k4, temp, '+');
    temp = vectorRDigit(Tchange / 6.0, temp, '*');
    temp = vectorOperation(y_i, temp, '+');
}

vector<TT> rungeKutta4Calc(const vector<TT> &yi, vector<TT> &tempVector, const TT tau, const TT Tchange) {
    vector<TT> resVector(yi.size());
    // calc k_i
    vector<TT> k1(f(yi));
    resVector = vectorOperation(tempVector, vectorRDigit(0.5 * Tchange * tau, k1, '*'), '+');
    vector<TT> k2(f(resVector));
    resVector = vectorOperation(tempVector, vectorRDigit(0.5 * Tchange * tau, k2, '*'), '+');
    vector<TT> k3(f(resVector));
    resVector = vectorOperation(tempVector, vectorRDigit(Tchange * tau, k3, '*'), '+');
    vector<TT> k4(f(resVector));

    // sum k_i
    k_iSum(resVector, k1, k2, k3, k4, yi, Tchange);

    return resVector;
}

vector<vector<TT>> rungeKutta2(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    vector<TT> k1(numOfPoints);
    vector<TT> k2(numOfPoints);
    for (int i = 0; i < n - 1; i++) {
        k1 = f(y[i]);
        vector<TT> temp = vectorOperation(y[i],
                                          vectorRDigit(step, k1, '*'), '+');
        k2 = f(temp);
        temp = vectorOperation(k1, k2, '+');
        y[i + 1] = vectorOperation(y[i], vectorRDigit((step / 2.0), temp, '*'), '+');
    }
    return y;
}

//vector<vector<TT>> rungeKutta4(const vector<TT> &cond, const TT Tchange, const TT eps = COMPARE_RATE) {
//    vector<vector<TT>> y(numOfPoints, vector<TT>(cond.size()));
//    y[0] = cond;
//
//    for (int i = 0; i < numOfPoints - 1; i++) {
//        vector<TT> defaultStepResult = rungeKutta4Calc(y[i], y[i], tau);
//
//        vector<TT> changedStepResult(defaultStepResult);
//        // auto step
//        for (int j = 0; j < numOfPoints - 1; j++) {
//            changedStepResult = rungeKutta4Calc(y[i], changedStepResult, tau / 2);
//        }
//        if (norm1Vector(vectorRDigit(1 / (pow(2, 4) - 1),
//                                     vectorOperation(changedStepResult, defaultStepResult, '-'), '*'))
//            > COMPARE_RATE) {
//            y[i + 1] = changedStepResult;
//        } else {
//            y[i + 1] = defaultStepResult;
//        }
//    }
//    return y;
//}

pair<vector<vector<TT>>, TT>
rungeKuttaRungeStep(const vector<TT> &cond, const TT Tchange, const TT eps = COMPARE_RATE) {
    vector<vector<TT>> y(numOfPoints, vector<TT>(cond.size()));
    y[0] = cond;
    TT tau = 1.0;
    const TT Tstep = Tchange / TT(numOfPoints);

    for (int i = 0; i < numOfPoints - 1; i++) {
        vector<TT> defaultStepResult = rungeKutta4Calc(y[i], y[i], tau, Tstep);

        vector<TT> changedStepResult(defaultStepResult);
        // auto step
        do {
            tau /= 2;
            defaultStepResult = changedStepResult;

            for (int j = 0; j < numOfPoints - 1; j++) {
                changedStepResult = rungeKutta4Calc(y[i], changedStepResult, tau, Tstep);
            }
        } while (norm1Vector(vectorRDigit(1 / (pow(2, 4) - 1),
                                          vectorOperation(changedStepResult, defaultStepResult, '-'), '*'))
                 > eps);

        y[i + 1] = changedStepResult;
        tau *= 2;
    }

    return {y, tau};
}

void templateOutput(const calcMethod method) {
    vector<vector<TT>> result;

    switch (method) {
        case MexplicitEuler:
            result = explicitEuler(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMexplicitEuler.txt");
            break;
        case MimplicitEuler:
            result = implicitEuler(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMimplicitEuler.txt");
            break;
        case Msymmetric:
            result = symmetric(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMsymmetric.txt");
            break;
        case MrungeKutta2:
            result = rungeKutta2(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMrungeKutta2.txt");
            break;
        case MrungeKutta4:
//            result = rungeKutta4(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMrungeKutta4.txt");
            break;
        case MrungeKuttaRungeStep:
            result = rungeKuttaRungeStep(initPoints, numOfPoints).first;
            outputMatrix(result, ADD_DOTS"data/outMrungeKuttaRungeStep.txt");
            break;
    }
    outputOnTheScreenMatrix(result);
    std::cout << "diff:\n";
    outputOnTheScreenMatrix(calcDiff(result));
}

void ImplicitEuler() {
    templateOutput(MexplicitEuler);
}

void ExplicitEuler() {
    templateOutput(MimplicitEuler);
}

void Symmetric() {
    templateOutput(Msymmetric);
}

void RungeKutta2() {
    templateOutput(MrungeKutta2);
}

void RungeKutta4() {
    templateOutput(MrungeKutta4);
}

void RungeKuttaChangeStep() {
    templateOutput(MrungeKuttaRungeStep);
}

TT calcDiff(const vector<vector<TT>> &answer, const TT Tstep) {
    vector<vector<TT>> diff(answer);//; = rungeKuttaRungeStep(initPoints, numOfPoints);

    for (int i = 0; i < answer.size(); ++i) {
        diff[i][0] = abs(cos(i * Tstep) - answer[i][0]);
        diff[i][1] = abs(-sin(i * Tstep) - answer[i][1]);
    }

//    return norm1Matrix(matrixOperations(ans2, answer, '-'));
    return norm1Matrix(diff);
}

void RungeKuttaGraph() {
    vector<vector<TT>> timeAndStep;
    timeAndStep.reserve(100);
    for (int i = 0; i < 100; ++i) {
        timeAndStep.push_back({TT((double) i / 10), rungeKuttaRungeStep(initPoints, TT((double) i / 10)).second});
    }
    outputMatrix(timeAndStep, ADD_DOTS"outTimeAndStep.txt");

    vector<vector<TT>> timeAndError;
    timeAndError.reserve(100);
    for (int i = 0; i < 100; ++i) {
        timeAndError.push_back(
                {TT((double) i / 10) / numOfPoints, calcDiff(rungeKuttaRungeStep(initPoints, TT((double) i / 10)).first,
                                                             (double) i / 10)});
    }
    outputMatrix(timeAndError, ADD_DOTS"outTimeAndError.txt");
}
