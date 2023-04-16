import { Epsilon, SmallEpsilon, SmallestEpsilon } from "../../../utils";
import Vector from "../../../vector";
import * as RootFinding from "./exports"


describe("Root finding: nonlinear", () => {
    test("1d", () => {

    });
    test("Nd", () => {
        let func = (p: Vector) => {
            const x = p.get(0);
            const y = p.get(1);
            return new Vector([20 * Math.log(x - y) - x - y - 6, 20 * Math.sin(0.7 * (x - y)) + 7 * (x + y)]);
        };
        let p0 = new Vector([0, -1]);
        let expectedRoot = new Vector([-0.46584782, -1.67846886]);
        let solver = new RootFinding.NewtonRaphson.Solver();
        let params = new RootFinding.NewtonRaphson.Params();
        params.fTolAbs = SmallEpsilon;
        params.fDotTolAbs = SmallEpsilon;
        params.jacobianEpsilon = SmallestEpsilon;
        let root: Vector;
        expect(root = RootFinding.NewtonRaphson.Solver.solve(func, p0, 10, params)).not.toThrow();
        expect(func(expectedRoot).lInfNorm()).toBeLessThan(params.fTolAbs);
        expect(Vector.sub(expectedRoot, root)).toBeLessThan(Epsilon);
    })
});

test.only("1D root finding", () => {
    const func = (x: number): number => { return x * x * x - 3.0 * x * x * Math.cos(x * 4) + Math.sin(3.3 * x); };

    let xa = -0.5;
    let xb = 1.0;
    let solution = 0.0;
    expect(RootFinding.Bisection.solve(func, xa, xb, 20)).toBeCloseTo(solution);
    expect(RootFinding.RegulaFalsi.solve(func, xa, xb, 20)).toBeCloseTo(solution);
});