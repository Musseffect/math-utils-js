
// todo: solvers for rect matrices
// todo: eigenvalue solvers

import Matrix from "../../denseMatrix";
import { Tolerance, SmallTolerance, assert, SmallestTolerance } from "../../utils";
import Vector from "../../vector";
import * as linSolvers from "./exports";

interface TestCase {
    m: Matrix;
    rhs: Vector;
    exactSolution: Vector;
    inverse: Matrix;
    determinant: number;
}

interface Tests {
    posDef: TestCase[],
    general: TestCase[]
};

const tests: Tests = { posDef: [], general: [] };

(function () {
    tests.posDef.push({
        m: new Matrix([10, -1, 2, 0, -1, 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8], 4, 4),
        rhs: new Vector([6, 25, -11, 15]),
        exactSolution: new Vector([1, 2, -1, 1]),
        inverse: new Matrix([259 / 2465, 23 / 2465, -3 / 145, -3 / 493, 23 / 2465, 758 / 7395, 2 / 435, -56 / 1479, -3 / 145, 2 / 435, 46 / 435, 1 / 87, -3 / 493, -56 / 1479, 1 / 87, 208 / 1479], 4, 4),
        determinant: 7395
    });
    /*tests.posDef.push({
        m: new Matrix([1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 3, 3, 1, 2, 3, 4], 4, 4),
        rhs: new Vector([-4, -7, -9, -10]),
        exactSolution: new Vector([-1, -1, -1, -1]),
        inverse: new Matrix([2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 1], 4, 4),
        determinant: 1
    });*/
    /*tests.posDef.push({
        m: new Matrix([4, 12, -16, 12, 37, -43, -16, -43, 98], 3, 3),
        rhs: new Vector([76, 215, -396]),
        exactSolution: new Vector([1, 2, -3]),
        inverse: new Matrix([1777 / 36, -122 / 9, 19 / 9, -122 / 9, 34 / 9, -5 / 9, 19 / 9, -5 / 9, 1 / 9], 3, 3),
        determinant: 36
    });*/
    tests.general.push({
        m: new Matrix([0.02, 0.01, 0, 0, 1, 2, 1, 0, 0, 1, 2, 1, 0, 0, 100, 200], 4, 4),
        rhs: new Vector([0.02, 1, 4, 800]),
        exactSolution: new Vector([1, 0, 0, 4]),
        inverse: new Matrix([80, -0.6, 0.4, -0.002, -60, 1.2, -0.8, 0.004, 40, -0.8, 1.2, -0.006, -20, 0.4, -0.6, 0.008], 4, 4),
        determinant: 5
    });
    /*
    tests.general.push({
        m: new Matrix([0, 1, 1, 1], 2, 2),
        rhs: new Vector([2, -1]),
        exactSolution: new Vector([-3, 2]),
        inverse: new Matrix([-1, 1, 1, 0], 2, 2),
        determinant: -1
    });*/
    const checkTest = (test: TestCase) => {
        assert(test.m.isSquare(), "Expected square matrix");
        assert(test.inverse.isSquare(), "Expected square inverse");
        assert(Vector.near(test.rhs, Matrix.postMulVec(test.m, test.exactSolution), SmallestTolerance), "Incorrect solution");
        assert(Vector.near(test.exactSolution, Matrix.postMulVec(test.inverse, test.rhs), SmallestTolerance), "Incorrect inverse");
    };
    for (const test of tests.posDef) {
        assert(test.m.isSymmetric(), "Expected symmetric matrix");
        assert(test.inverse.isSymmetric(), "Expected symmetric inverse");
        checkTest(test);
    }
    for (const test of tests.general) {
        checkTest(test);
    }
})();

describe('Linear solvers (dense square matrices)', () => {
    describe.skip.each(tests.posDef)('Symmetric positive definite matrices %#', (testCase: TestCase) => {
        expect(testCase.m.isSymmetric()).toBeTruthy();
        describe('Factorizations', () => {
            test('LL', () => {
                let solver = new linSolvers.LLT(null, SmallTolerance);
                expect(() => solver.factorize(testCase.m)).not.toThrow();
                expect(solver.LLT).not.toBeNull();
                console.log(`m: ${testCase.m.toString()}`);
                console.log(`LLT: ${solver.LLT.toString()}`);
                expect(Vector.sub(solver.solve(testCase.rhs) as Vector, testCase.exactSolution).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
                expect(Matrix.lInfDistance(testCase.inverse, solver.inverse() as Matrix)).toBeLessThan(SmallTolerance);
                expect(solver.determinant()).toBeCloseTo(testCase.determinant, 4);
            });
            test.skip('LDL', () => {
                let solver = new linSolvers.LDLT(null, SmallTolerance);
                expect(() => solver.factorize(testCase.m)).not.toThrow();
                expect(solver.LDLT).not.toBeNull();
                console.log(`m: ${testCase.m.toString()}`);
                console.log(`LDLT: ${solver.LDLT.toMatrix().toString()}`);
                expect(Vector.sub(solver.solve(testCase.rhs) as Vector, testCase.exactSolution).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
                expect(Matrix.lInfDistance(testCase.inverse, solver.inverse() as Matrix)).toBeLessThan(SmallTolerance);
                expect(solver.determinant()).toBeCloseTo(testCase.determinant, 4);
            });
            test('PPLU', () => {
                let solver = new linSolvers.PartialPivLU(null, SmallTolerance);
                expect(() => solver.factorize(testCase.m)).not.toThrow();
                expect(solver.LU).not.toBeNull();
                expect(Vector.sub(solver.solve(testCase.rhs) as Vector, testCase.exactSolution).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
                expect(Matrix.lInfDistance(testCase.inverse, solver.inverse() as Matrix)).toBeLessThan(SmallTolerance);
                expect(solver.determinant()).toBeCloseTo(testCase.determinant, 4);
            });
            test.skip('FPLU', () => {
                let solver = new linSolvers.FullPivLU(null, SmallTolerance);
                expect(() => solver.factorize(testCase.m)).not.toThrow();
                expect(solver.LU).not.toBeNull();
                expect(Vector.sub(solver.solve(testCase.rhs) as Vector, testCase.exactSolution).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
                expect(Matrix.lInfDistance(testCase.inverse, solver.inverse() as Matrix)).toBeLessThan(SmallTolerance);
                expect(solver.determinant()).toBeCloseTo(testCase.determinant, 4);
            });
        });
        describe.skip('Iterative', () => {
            test('ConjGrad', () => {
                let result = linSolvers.CG.solve(testCase.m, testCase.rhs);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
            test('gauss-zeidel', () => {
                let result = linSolvers.gaussSeidel.solve(testCase.m, testCase.rhs, 30, SmallTolerance);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
            test('jacobi', () => {
                let result = linSolvers.jacobi.solve(testCase.m, testCase.rhs, 30, SmallTolerance);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
            test('sor', () => {
                let result = linSolvers.sor.solve(testCase.m, testCase.rhs, 30, 1.0, SmallTolerance);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
        });
    });
    describe.each(tests.general)('General matrices', (testCase: TestCase) => {
        describe('Factorizations', () => {
            test.skip('PartialPivLU', () => {
                /*expect(Vector.near(linSolvers.PartialPivLU.solve(testCase.m.clone(), testCase.rhs.clone(), SmallTolerance), testCase.exactSolution, Tolerance)).toBeTruthy();
                // test decomposition
                let luSolver = new linSolvers.PartialPivLU(testCase.m);
                let PT = luSolver.P.inverse().toMatrix();
                let PTLU = Matrix.mul(PT, Matrix.mul(luSolver.L.toMatrix(), luSolver.U.toMatrix()));
                expect(Math.abs(PT.determinantNaive())).toBeCloseTo(1);
                expect(Matrix.near(PTLU, testCase.m));

                expect(Matrix.near(linSolvers.PartialPivLU.solveMatrix(testCase.m.clone(), Matrix.identity(testCase.m.height())), testCase.m.inverseNaive()));
*/
                let solver = new linSolvers.PartialPivLU(null, SmallTolerance);
                expect(() => solver.factorize(testCase.m)).not.toThrow();
                expect(solver.LU).not.toBeNull();
                expect(Vector.sub(solver.solve(testCase.rhs) as Vector, testCase.exactSolution).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
                expect(Matrix.lInfDistance(testCase.inverse, solver.inverse() as Matrix)).toBeLessThan(SmallTolerance);
                expect(solver.determinant()).toBeCloseTo(testCase.determinant, 4);
            });
            test('FullPivLU', () => {
                /*// test factorization
                let lusolver = new linSolvers.FullPivLU(testCase.m);
                let PT = luSolver.P.inverse().toMatrix();
                let QT = luSolver.Q.inverse().toMatrix();
                let PTLUQT = Matrix.mul(Matrix.mul(PT, Matrix.mul(luSolver.L.toMatrix(), luSolver.U.toMatrix())), QT);
                expect(Math.abs(PT.determinantNaive())).toBeCloseTo(1);
                expect(Math.abs(QT.determinantNaive())).toBeCloseTo(1);
                expect(Matrix.near(PTLUQT, testCase.m));
            
                expect(Matrix.near(linSolvers.PartialPivLU.solveMatrix(testCase.m.clone(), Matrix.identity(testCase.m.height())), testCase.m.inverseNaive()));
            */
                let solver = new linSolvers.FullPivLU(null, SmallTolerance);
                expect(() => solver.factorize(testCase.m)).not.toThrow();
                console.log(`p ${solver.P.array()}`);
                console.log(`q ${solver.Q.array()}`);
                expect(solver.LU).not.toBeNull();
                expect(Vector.sub(solver.solve(testCase.rhs) as Vector, testCase.exactSolution).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
                expect(Matrix.lInfDistance(testCase.inverse, solver.inverse() as Matrix)).toBeLessThan(SmallTolerance);
                expect(solver.determinant()).toBeCloseTo(testCase.determinant, 4);
            });

        });
        describe.skip('Iterative', () => {
            test('gauss-zeidel', () => {
                let result = linSolvers.gaussSeidel.solve(testCase.m, testCase.rhs, 30, SmallTolerance);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
            test('jacobi', () => {
                let result = linSolvers.jacobi.solve(testCase.m, testCase.rhs, 30, SmallTolerance);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
            test('sor', () => {
                let result = linSolvers.sor.solve(testCase.m, testCase.rhs, 30, 1.0, SmallTolerance);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
        });
    });
});

test.skip('Linear solvers (dense)', () => {
    interface testData {
        m: Matrix;
        rhs: Vector;
        exactSolution: Vector;
    };
    let testExamples: testData[] = [];
    let m = Matrix.empty(4, 4);
    let rhs = new Vector([6, 25, -11, 15]);
    m.set(0, 0, 10);
    m.set(0, 1, -1);
    m.set(0, 2, 2);

    m.set(1, 0, -1);
    m.set(1, 1, 11);
    m.set(1, 2, -1);
    m.set(1, 3, 3);

    m.set(2, 0, 2);
    m.set(2, 1, -1);
    m.set(2, 2, 10);
    m.set(2, 3, -1);

    m.set(3, 1, 3);
    m.set(3, 2, -1);
    m.set(3, 3, 8);
    let exactSolution = new Vector([1, 2, -1, 1]);
    testExamples.push({ m, rhs, exactSolution });
    for (let testExample of testExamples) {
        console.log(`Matrix ${testExample.m.toString()}`);
        console.log(`Determinant ${testExample.m.determinantNaive()}`);
        const inverse = testExample.m.inverseNaive();
        expect(Vector.near(linSolvers.gaussSeidel.solve(testExample.m, testExample.rhs, 15, SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();
        expect(Vector.near(linSolvers.jacobi.solve(testExample.m, testExample.rhs, 25, SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();
        expect(Vector.near(linSolvers.sor.solve(testExample.m, testExample.rhs, 35, 0.5, SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();
        expect(Vector.near(linSolvers.PartialPivLU.solve(testExample.m, testExample.rhs, SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();
        expect(Vector.near(linSolvers.FullPivLU.solve(testExample.m, testExample.rhs, SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();

        let identity = Matrix.identity(testExample.m.width());
        expect(Matrix.near(Matrix.mul(testExample.m, inverse), identity)).toBeTruthy();
        expect(Matrix.near(linSolvers.PartialPivLU.solveMatrix(testExample.m, identity, SmallTolerance), inverse, Tolerance)).toBeTruthy();
        console.log(`Inverse ${inverse.toString()}`);
        expect(Matrix.near(linSolvers.FullPivLU.solveMatrix(testExample.m, identity, SmallTolerance), inverse, Tolerance)).toBeTruthy();

        const testDimensions = (solution: Matrix, A: Matrix, Rhs: Matrix) => {
            expect(solution.width() == Rhs.width());
            expect(solution.height() == A.width());
        };
        const Rhs = Matrix.empty(testExample.m.height(), 2);
        testDimensions(
            linSolvers.PartialPivLU.solveMatrix(testExample.m, Rhs), testExample.m, Rhs);
        testDimensions(
            linSolvers.FullPivLU.solveMatrix(testExample.m, Rhs), testExample.m, Rhs);
    }


    m = Matrix.empty(3, 3);
    m.set(0, 0, 4);
    m.set(0, 1, 12);
    m.set(0, 2, -16);

    m.set(1, 0, 12);
    m.set(1, 1, 37);
    m.set(1, 2, -43);

    m.set(2, 0, -16);
    m.set(2, 1, -43);
    m.set(2, 2, 98);
    exactSolution = new Vector([1, 2, -3]);
    rhs = Matrix.postMulVec(m, exactSolution);
    expect(Vector.near(linSolvers.LLT.solve(m.clone(), rhs), exactSolution, Tolerance)).toBeTruthy();

    /*
        let singularMatrix = Matrix.empty(3, 3);
        // rank 2 singular matrix
        singularMatrix.set(1, 1, 1);
        singularMatrix.set(2, 2, 1);
    */
});

test.skip("QR tests", () => {


});