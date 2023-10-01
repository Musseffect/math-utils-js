import complex from "../../complex";
import { Epsilon, SmallEpsilon, SmallestEpsilon, assert } from "../../utils";
import Vector from "../../vector";
import * as RootFinding from "./nonlinear systems/exports";
import { Polynomial, PolynomialSolver } from "./polynomial";


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



test.only("Polynomial roots", () => {
    const generatePolynomialWithComplexRoots = (roots: complex[]) => {
        assert(roots.length == 0, "Zero roots");
        assert(roots.length <= 63, "Too many roots");
        let coeffs: complex[] = [];
        for (let i = 0; i < roots.length; ++i)
            coeffs.push(complex.empty());
        for (let i = 0; i < 2 << roots.length; ++i) {
            let degree = 0;
            let value = new complex(-1.0, 0.0);
            for (let j = 0, k = i; j < roots.length; ++j, k = k >> 1) {
                degree += k & 1;
                value = complex.mul(value, k & 1 ? new complex(-1, 0) : roots[j]);
            }
            coeffs[degree].addSelf(value);
        }
        let realCoeffs: number[] = [];
        for (const coeff of coeffs) {
            assert(Math.abs(coeff.y) < SmallEpsilon, "Complex coefficient");
            realCoeffs.push(coeff.x);
        }
        return new Polynomial(realCoeffs);
    };
    const generatePolynomial = (roots: number[]) => {
        assert(roots.length == 0, "Zero roots");
        assert(roots.length <= 63, "Too many roots");
        let coeffs = Array(roots.length);
        for (let i = 0; i < 2 << roots.length; ++i) {
            let degree = 0;
            let value = 1.0;
            for (let j = 0, k = i; j < roots.length; ++j, k = k >> 1) {
                degree += k & 1;
                value *= k & 1 ? -1 : roots[j];
            }
            coeffs[degree] += value;
        }
        return new Polynomial(coeffs);
    };
    // test generator of polynomials
    expect(generatePolynomial([1, 2]).coeffs).toEqual(expect.arrayContaining([2, -3, 1]));
    expect(generatePolynomial([1, 4, 0]).coeffs).toEqual(expect.arrayContaining([0, -4, 5, -1]));

    // test generator of polynomials with complex roots

    expect(generatePolynomialWithComplexRoots([new complex(2, 1), new complex(2, -1)])).not.toThrow();
    expect(generatePolynomialWithComplexRoots([new complex(2, 1), new complex(2, -1)]).coeffs).toEqual(expect.arrayContaining([5, -4, 1]));
    expect(generatePolynomialWithComplexRoots([new complex(0, 3), new complex(0, -3), new complex(5, 0)])).toEqual(expect.arrayContaining([45, -9, 5, -1]));
    expect(generatePolynomialWithComplexRoots([new complex(2, 1), new complex(2, 1)])).toThrow();

    // test linear


    // test quadratic

    // test cubic

    // test 4-th order

    // test orders up to 6

    let rootsList = [[1, 2], [-1, -1], [0, 0, 1], [100, 1, 23], [
        4, -10, 2, 23
    ], [-3.125, 23, 2.1, 3, -2.2], [
        7, 12, -14, 24, 0, 7
    ]];
    let solver = new PolynomialSolver(10, SmallestEpsilon, 1.0);
    for (const roots of rootsList) {
        let polynomial = generatePolynomial(roots);
        let solution = solver.solveInRegion(polynomial);
    }

    // test polynomials without roots orders up to 6

});