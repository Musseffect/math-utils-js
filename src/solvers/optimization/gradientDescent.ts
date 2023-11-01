import { SmallTolerance } from "../../utils";
import Vector from "../../vector";
import { LineSearch } from "../line search/lineSearch";
import { OptimizationProblem } from "./optimizationProblem";

enum LineSearchAlgorithm {
    Wolf = 0,
    Armijo = 1,
    Goldstein = 2
    // todo: add more
};

class GradientDescent {
    lineSearch: LineSearch;
    tolerance: number;
    step: number;
    constructor(lineSearchAlg: LineSearchAlgorithm = LineSearchAlgorithm.Wolf, tolerance: number = SmallTolerance) {
        //this.lineSearch = lineSearchAlg;
        this.tolerance = tolerance;
    };
    solve(f: OptimizationProblem, x0: Vector, numIters: number): Vector {
        let x = x0.clone();
        for (let i = 0; i < numIters; ++i) {
            let direction = Vector.negate(f.dfdx(x));
            if (direction.l2Norm() < this.tolerance) return x;
            let step = this.lineSearch.step(x, direction, this.step);
            x.addSelf(direction.scaleSelf(step));
        }
        throw new Error("Not implemented");
    }
}