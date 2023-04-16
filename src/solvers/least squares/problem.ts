import Vector from "../../vector";

export abstract class LeastSquaresFunction {
    abstract numParameters(): number;
    abstract numVariables(): number;
    abstract f(x: Vector, params: Vector): number;
}
export class LeastSquaresProblem {
    xPoints: Vector[];
    yValues: number[];
}