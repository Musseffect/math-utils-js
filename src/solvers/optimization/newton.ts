import { SmallTolerance } from "../../utils";
import Vector from "../../vector";
import LDLT from "../linear systems/ldlt";
import { OptimizationProblem } from "./optimizationProblem";


export class Newton {
    static solve(op: OptimizationProblem, x0: Vector, numIters: number, xAbsTol: number = SmallTolerance, dfAbsTol: number = SmallTolerance) {
        let x = x0.clone();
        let fValue = op.f(x);
        for (let i = 0; i < numIters; ++i) {
            let grad = op.dfdx(x);
            if (grad.l2Norm() > dfAbsTol) return x;
            let hessian = op.dfdxdy(x);
            // todo: regularize
            let ldlt = new LDLT(hessian);
        }
    }
    // calc H, calc QLQ decomposition and clamp negative eigenvalues or try some other scheme
}