
import { complex } from "./complex";
import { Polynomial, PolynomialSolver, generatePolynomial, generatePolynomialWithComplexRoots } from "./polynomial";
import { calcEigenvalues } from "./solvers/linear systems/eigenvalues";
import { SmallTolerance, SmallestTolerance, near } from "./utils";




test.skip("Polynomial operations", () => {
    // addition
    // multiplication
    // division
    // derivative
    // integration
});


test("Polynomial roots", () => {

    // test generator of polynomials
    expect(generatePolynomial([1, 2]).coeffs).toEqual(expect.arrayContaining([2, -3, 1]));
    expect(generatePolynomial([1, 3, 5]).coeffs).toEqual(expect.arrayContaining([15, -23, 9, -1]));
    expect(generatePolynomial([1, 4, 5, 0]).coeffs).toEqual(expect.arrayContaining([0, -20, 29, -10, 1]));

    // test generator of polynomials with complex roots

    expect(() => generatePolynomialWithComplexRoots([new complex(2, 1), new complex(2, -1)])).not.toThrow();
    expect(generatePolynomialWithComplexRoots([new complex(2, 1), new complex(2, -1)]).coeffs).toEqual(expect.arrayContaining([5, -4, 1]));
    expect(generatePolynomialWithComplexRoots([new complex(0, 3), new complex(0, -3), new complex(5, 0)]).coeffs).toEqual(expect.arrayContaining([45, -9, 5, -1]));
    expect(() => generatePolynomialWithComplexRoots([new complex(2, 1), new complex(2, 1)])).toThrow();

    // test linear


    // test quadratic

    // test cubic

    // test 4-th order

    // test orders up to 6

    let rootsList = [/*[1, 2], [-1, -1], [0, 0, 1], [100, 1, 23], [10, 10, 3], [
        4, -10, 2, 23
    ], [-3.125, 23, 2.1, 3, -2.2],*/ [
            7, 12, -14, 24, 0, 7
        ]];
    let solver = new PolynomialSolver(50, SmallestTolerance, 0.5);
    for (const roots of rootsList) {
        let polynomial = generatePolynomial(roots);
        let solution = solver.solveInRegion(polynomial);
        if (solution.length != roots.length) {
            console.log(`Polynomial: ${polynomial.toString()}`);
            let der = polynomial;
            for (let i = 0; i < polynomial.degree() - 1; ++i) {
                der = der.derivative();
                console.log(`Der ${i}: ${der.toString()}`);
            }
        }
        expect(solution.length).toBeCloseTo(roots.length);
        let sortedRoots = roots.sort((a, b) => a - b);
        let error: number[] = [];
        console.log(`Expected: ${sortedRoots.toString()}, actual: ${solution.toString()}`);
        for (let i = 0; i < roots.length; ++i) {
            expect(Math.abs(sortedRoots[i] - solution[i])).toBeLessThanOrEqual(SmallTolerance);
            error.push(Math.abs(sortedRoots[i] - solution[i]));
        }
        console.log(`error: ${error.toString()}`);
    }

    // test polynomials without roots orders up to 6


    // generate matrix and solve eigenvalues
    expect(calcEigenvalues(generatePolynomial([1, 4, 0]).companionMatrix(), 10, SmallestTolerance)).toEqual(expect.arrayContaining([1, 4, 5, 0]));
});