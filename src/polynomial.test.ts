
import { complex } from "./complex";
import { Polynomial, PolynomialSolver, generatePolynomial, generatePolynomialWithComplexRoots } from "./polynomial";
import { calcEigenvalues } from "./solvers/linear systems/eigenvalues";
import { SmallestTolerance } from "./utils";




test.only("Polynomial roots", () => {

    // test generator of polynomials
    expect(generatePolynomial([1, 2]).coeffs).toEqual(expect.arrayContaining([2, -3, 1]));
    expect(generatePolynomial([1, 4, 5, 0]).coeffs).toEqual(expect.arrayContaining([0, -4, 5, -1]));

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
    let solver = new PolynomialSolver(10, SmallestTolerance, 1.0);
    for (const roots of rootsList) {
        let polynomial = generatePolynomial(roots);
        let solution = solver.solveInRegion(polynomial);
    }

    // test polynomials without roots orders up to 6


    // generate matrix and solve eigenvalues
    expect(calcEigenvalues(generatePolynomial([1, 4, 0]).companionMatrix(), 10, SmallestTolerance)).toEqual(expect.arrayContaining([1, 4, 5, 0]));
});