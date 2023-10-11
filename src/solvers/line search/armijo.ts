import { assert, Tolerance, SmallTolerance, SmallestTolerance } from "../../utils";
import Vector from "../../vector";
import { LineSearch, LineSearchProblem } from "./lineSearch";

export class ArmijoBacktracking extends LineSearch {
    tau: number = 0.5;
    c: number = 0.5;
    numIters = 100;
    tolerance = SmallTolerance;
    constructor(problem: LineSearchProblem) {
        super(problem);
    }
    setTau(value: number) {
        assert(value > 0 && value < 1.0, "Invalid value for parameter tau");
        this.tau = value;
    }
    setC(value: number) {
        assert(value > 0 && value < 1.0, "Invalid value for parameter c");
        this.c = value;
    }
    stMaxIters(numIters: number) {
        this.numIters = numIters;
    }
    setTolerance(tolerance: number) {
        this.tolerance = tolerance;
    }
    public step(x: Vector, direction: Vector, initialStep: number = 1.0): number {
        let initialValue = this.problem.f(x);
        let step = initialStep;
        const m = Vector.dot(this.problem.grad(x), direction);
        if (m > -SmallestTolerance) return this.tolerance;
        const t = -this.c * m;
        for (let iter = 0; iter < this.numIters && step > this.tolerance; ++iter) {
            let curValue = this.problem.f(Vector.add(x, Vector.scale(direction, step)));
            if ((initialValue - curValue) >= step * t) break;
            step = this.tau * step;
        }
        return step;
    }
}

class ArmijoTwoWayBacktracing extends ArmijoBacktracking {
    largestStep: number;
    setLargestStep(value: number) {
        this.largestStep = value;
    }
    public step(x: Vector, direction: Vector, initialStep: number = 1.0): number {
        let initialValue = this.problem.f(x);
        let step = initialStep;
        const m = Vector.dot(this.problem.grad(x), direction);
        if (m > -SmallestTolerance) return this.tolerance;
        const t = -this.c * m;
        assert(t < 0, "Invalid direction");
        let iter = 0;
        for (; iter < this.numIters && step <= this.largestStep; ++iter) {
            let curValue = this.problem.f(Vector.add(x, Vector.scale(direction, step)));
            if (initialValue - curValue < step * t) break;
            step = step / this.tau;
        }
        for (; iter < this.numIters && step >= this.tolerance; ++iter) {
            let curValue = this.problem.f(Vector.add(x, Vector.scale(direction, step)));
            if (initialValue - curValue >= step * t) break;
            step = step * this.tau;
        }
        return step;
    }
}
