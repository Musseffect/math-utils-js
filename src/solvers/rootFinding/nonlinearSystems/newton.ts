import matrix from "../../../denseMatrix";
import { SmallEpsilon } from "../../../utils";
import vector from "../../../vector";
import { ConvergenseFailureException } from "../linearSystems/exceptions";
import FullPivLU from "../linearSystems/fullPivLU";


const SolverName = "Newton";
class Params {
    jacobian?: (x: vector) => matrix;
    fDotTolAbs: number = SmallEpsilon;
    fTolAbs: number = SmallEpsilon;
    jacobianEpsilon: number = SmallEpsilon;
}

class Solver {
    private static numericJacobian(f: (x: vector) => vector, x: vector, step: number): matrix {
        throw Error("Not implemented");
    }
    static solve(f: (x: vector) => vector, x0: vector, numIters: number, params:Params = new Params()):vector {
        let x = x0.clone();
        let fVec = vector.negate(f(x));
        let maxNorm = fVec.lInfNorm();
        for (let iter = 0; iter < numIters; ++iter) {
            if (maxNorm < params.fTolAbs) return x;
            let J: matrix;
            if (params.jacobian)
                J = params.jacobian(x);
            else
                J = numericJacobian(f, x, SmallEpsilon);
            if (J.maxNorm() < params.fDotTolAbs)
                return x;
            try {
                let dx = FullPivLU.solve(J, fVec);
                x.addSelf(dx);
            } catch (exc) {
                break;
            }
            fVec = vector.negate(f(x));
            maxNorm = fVec.lInfNorm();
        }
        throw new ConvergenseFailureException(SolverName);
    }
    // add line search, bounds, etc.
    // trust region methods
}

export default {Solver, Params};