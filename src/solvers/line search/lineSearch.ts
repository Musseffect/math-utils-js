import Matrix from "../../denseMatrix";
import Vector from "../../vector";

export abstract class LineSearchProblem {
    public abstract f(x: Vector): number;
    public abstract grad(x: Vector): Vector;
    public abstract hessian(x: Vector): Matrix;
}

export abstract class LineSearch {
    problem: LineSearchProblem;
    constructor(problem: LineSearchProblem) {
        this.problem = problem;
    }
    public abstract step(x: Vector, direction: Vector, initialStep: number): number;
}