import { SmallTolerance } from "../../utils";
import Vector from "../../vector";
import { ConvergenseFailureException } from "../linear systems/exceptions";
import LDLT from "../linear systems/ldlt";
import PartialPivLU from "../linear systems/partialPivLU";
import { Bounds, OptimizationProblem } from "./optimizationProblem";
import { LineSearchAlgorithm, initializeLineSearch } from "./utils";


export class Newton {
    lineSearchAlgo: LineSearchAlgorithm;
    initialStep: number;
    iterations: number;
    gradTolerance: number = SmallTolerance;
    dxAbsTol: number = SmallTolerance;
    constructor() {

    }
    solve(op: OptimizationProblem, x0: Vector, bounds?: Bounds): Vector {
        const lineSearch = initializeLineSearch(this.lineSearchAlgo, op);
        let x = x0.clone();
        let iter = 0;
        while (true) {
            let grad = op.dfdx(x);
            if (grad.l2Norm() < this.gradTolerance) return x;
            if (++iter == this.iterations) break;
            let hessian = op.dfdxdy(x);
            // todo: regularize
            let ldlt = new LDLT(hessian);
            throw new Error("Not implemented");
            let direction: Vector = PartialPivLU.solve(hessian, Vector.negate(grad));
            let step = this.initialStep;
            if (!bounds) {
                step = Math.min(step, bounds.intersect(x, direction));
            }
            step = lineSearch.step(x, direction, step);
            x.addSelf(direction.scaleSelf(step));
        }
        throw new ConvergenseFailureException("Newton");
    }
    // calc H, calc QLQ decomposition and clamp negative eigenvalues or try some other scheme
}