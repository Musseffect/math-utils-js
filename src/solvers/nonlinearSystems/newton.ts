import matrix from "../../denseMatrix";
import { SmallEpsilon } from "../../utils";
import vector from "../../vector";
import { ConvergenseFailureException } from "../linearSystems/exceptions";
import FullPivLU from "../linearSystems/fullPivLU";


const Solver = "Newton";

export default class Newton {
    private static numericJacobian(f: (x: vector) => vector, x: vector, step: number): matrix {

    }
    static solve(f: (x: vector) => vector, x0: vector, numIters: number, fTolAbs: number, jacobian?: (x: vector) => matrix): vector {
        let x = x0.clone();
        let fVec = vector.negate(f(x));
        let maxNorm = fVec.lInfNorm();
        for (let iter = 0; iter < numIters; ++iter) {
            if (maxNorm < fTolAbs) return x;
            let J: matrix;
            if (jacobian)
                J = jacobian(x);
            else
                J = numericJacobian(f, x, SmallEpsilon);
            try {
                let dx = FullPivLU.solve(J, fVec);
                x.addSelf(dx);
            } catch (exc) {
                break;
            }
            fVec = vector.negate(f(x));
            maxNorm = fVec.lInfNorm();
        }
        throw new ConvergenseFailureException(Solver);
    }
}