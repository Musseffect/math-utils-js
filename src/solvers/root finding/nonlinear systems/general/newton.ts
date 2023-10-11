import Matrix from "../../../../denseMatrix";
import { assert, SmallTolerance } from "../../../../utils";
import vector from "../../../../vector";
import { ConvergenseFailureException } from "../../../linear systems/exceptions";
import FullPivLU from "../../../linear systems/fullPivLU";


const SolverName = "Newton";
class Params {
    jacobian?: (x: vector) => Matrix;
    fDotTolAbs: number = SmallTolerance;
    fTolAbs: number = SmallTolerance;
    jacobianEpsilon: number = SmallTolerance;
}

class Solver {
    private static numericJacobian(f: (x: vector) => vector, x: vector, step: number): Matrix {
        let result: Matrix = Matrix.empty(x.size(), x.size());
        // forward difference
        let f0 = f(x);
        for (let i = 0; i < x.size(); ++i) {
            let scalar = x.get(i);
            x.set(i, scalar + step);
            let dfdx = f(x).subSelf(f0).scaleSelf(1.0 / step);
            x.set(i, scalar);
            result.setColumn(i, dfdx);
        }
        return result;
    }
    static solve(f: (x: vector) => vector, x0: vector, numIters: number, params: Params = new Params()): vector {
        let x = x0.clone();
        let fVec = vector.negate(f(x));
        let maxNorm = fVec.lInfNorm();
        let calcJ = params.jacobian ? params.jacobian : (y: vector) => { return Solver.numericJacobian(f, y, SmallTolerance) };
        for (let iter = 0; iter < numIters; ++iter) {
            if (maxNorm < params.fTolAbs) return x;
            let J: Matrix = calcJ(x);
            if (J.lInfNorm() < params.fDotTolAbs)
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
    // todo: add line search, bounds, etc.
    // trust region methods
    // check rank
}

export default { Solver, Params };