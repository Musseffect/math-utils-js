import AbstractMatrix from "../../../abstractMatrix";
import { SparseMatrixCSR } from "../../../sparseMatrix";
import Triplet from "../../../triplet";
import { SmallTolerance, SmallestTolerance, Tolerance, assert } from "../../../utils";
import Vector from "../../../vector";
import gaussSeidel from "./gaussSeidel";
import jacobi from "./jacobi";
import { ConjugateGradients } from "./conjugateGradients";
import Matrix from "../../../denseMatrix";

interface TestCase {
    m: SparseMatrixCSR;
    rhs: Vector;
    exactSolution: Vector;
    determinant: number;
    inverse: Matrix;
}

interface Tests {
    posDef: TestCase[],
    general: TestCase[]
};

const tests: Tests = { posDef: [], general: [] };

(function () {
    let triplets: Triplet[] = [];
    triplets.push({ row: 0, column: 0, value: 4 });
    triplets.push({ row: 0, column: 1, value: 1 });
    triplets.push({ row: 0, column: 2, value: 2 });
    triplets.push({ row: 0, column: 3, value: 0.5 });
    triplets.push({ row: 0, column: 4, value: 2 });
    triplets.push({ row: 1, column: 0, value: 1 });
    triplets.push({ row: 1, column: 1, value: 0.5 });
    triplets.push({ row: 2, column: 0, value: 2 });
    triplets.push({ row: 2, column: 2, value: 3 });
    triplets.push({ row: 3, column: 0, value: 0.5 });
    triplets.push({ row: 3, column: 3, value: 0.625 });
    triplets.push({ row: 4, column: 0, value: 2 });
    triplets.push({ row: 4, column: 4, value: 16 });
    tests.posDef.push({
        m: SparseMatrixCSR.fromTriplets(triplets, 5, 5),
        rhs: new Vector([6, 25, -11, 15, 10]),
        exactSolution: new Vector([-2995, 6040, 1993, 2420, 375]),
        determinant: 0.5,
        inverse: new Matrix([
            60, -120, -40, -48, -7.5,
            -120, 242, 80, 96, 15,
            -40, 80, 27, 32, 5,
            - 48, 96, 32, 40, 6,
            -7.5, 15, 5, 6, 1
        ], 5, 5)
    });
    triplets = [];
    triplets.push({ row: 0, column: 0, value: 2 });
    triplets.push({ row: 0, column: 1, value: -1 });
    triplets.push({ row: 1, column: 0, value: -1 });
    triplets.push({ row: 1, column: 1, value: 2 });
    triplets.push({ row: 1, column: 2, value: -1 });
    triplets.push({ row: 2, column: 1, value: -1 });
    triplets.push({ row: 2, column: 2, value: 2 });
    triplets.push({ row: 2, column: 3, value: -1 });
    triplets.push({ row: 3, column: 2, value: -1 });
    triplets.push({ row: 3, column: 3, value: 2 });
    triplets.push({ row: 3, column: 4, value: -1 });
    triplets.push({ row: 4, column: 3, value: -1 });
    triplets.push({ row: 4, column: 4, value: 2 });
    triplets.push({ row: 4, column: 5, value: -1 });
    triplets.push({ row: 5, column: 4, value: -0.5 });
    triplets.push({ row: 5, column: 5, value: 1 });
    triplets.push({ row: 5, column: 6, value: -0.5 });
    triplets.push({ row: 6, column: 5, value: -1 });
    triplets.push({ row: 6, column: 6, value: 2 });
    triplets.push({ row: 6, column: 7, value: -1 });
    triplets.push({ row: 7, column: 6, value: -0.5 });
    triplets.push({ row: 7, column: 7, value: 1 });
    triplets.push({ row: 7, column: 8, value: -0.5 });
    triplets.push({ row: 8, column: 6, value: 0.5 });
    triplets.push({ row: 8, column: 7, value: -1.5 });
    triplets.push({ row: 8, column: 8, value: 4 });
    triplets.push({ row: 8, column: 9, value: -1.5 });
    triplets.push({ row: 9, column: 8, value: -1 });
    triplets.push({ row: 9, column: 9, value: 2 });
    tests.general.push({
        m: SparseMatrixCSR.fromTriplets(triplets, 10, 10),
        rhs: new Vector([1, 3, 4, -4, 2, 3, 1, -2, -3, 1]),
        exactSolution: new Vector([
            514 / 83,
            945 / 83,
            1127 / 83,
            977 / 83,
            1159 / 83,
            1175 / 83,
            693 / 83,
            128 / 83,
            -105 / 83,
            -11 / 83
        ]),
        determinant: 10.375,
        inverse: new Matrix([
            74 / 83, 65 / 83, 56 / 83, 47 / 83, 38 / 83, 58 / 83, 20 / 83, 26 / 83, 4 / 83, 3 / 83,
            65 / 83, 130 / 83, 112 / 83, 94 / 83, 76 / 83, 116 / 83, 40 / 83, 52 / 83, 8 / 83, 6 / 83,
            56 / 83, 112 / 83, 168 / 83, 141 / 83, 114 / 83, 174 / 83, 60 / 83, 78 / 83, 12 / 83, 9 / 83,
            47 / 83, 94 / 83, 141 / 83, 188 / 83, 152 / 83, 232 / 83, 80 / 83, 104 / 83, 16 / 83, 12 / 83,
            38 / 83, 76 / 83, 114 / 83, 152 / 83, 190 / 83, 290 / 83, 100 / 83, 130 / 83, 20 / 83, 15 / 83,
            29 / 83, 58 / 83, 87 / 83, 116 / 83, 145 / 83, 348 / 83, 120 / 83, 156 / 83, 24 / 83, 18 / 83,
            20 / 83, 40 / 83, 60 / 83, 80 / 83, 100 / 83, 240 / 83, 140 / 83, 182 / 83, 28 / 83, 21 / 83,
            11 / 83, 22 / 83, 33 / 83, 44 / 83, 55 / 83, 132 / 83, 77 / 83, 208 / 83, 32 / 83, 24 / 83,
            2 / 83, 4 / 83, 6 / 83, 8 / 83, 10 / 83, 24 / 83, 14 / 83, 68 / 83, 36 / 83, 27 / 83,
            1 / 83, 2 / 83, 3 / 83, 4 / 83, 5 / 83, 12 / 83, 7 / 83, 34 / 83, 18 / 83, 55 / 83
        ], 10, 10)
    });
    const checkTest = (test: TestCase, testID: number) => {
        assert(test.m.isSquare(), `Expected square matrix ${testID}`);
        assert(test.exactSolution.size() == test.m.numCols(), `Inconsistent system size ${testID}`);
        assert(test.rhs.size() == test.m.numCols(), `Inconsistent system size ${testID}`);
        assert(Vector.near(test.rhs, SparseMatrixCSR.postMul(test.m, test.exactSolution), SmallestTolerance), `Incorrect solution ${testID}`);
        assert(Vector.near(test.exactSolution, Matrix.postMulVec(test.inverse, test.rhs), SmallestTolerance), `Incorrect inverse ${testID}`);
    };
    let testID = 0;
    for (const test of tests.posDef) {
        assert(test.m.isSymmetric(), "Expected symmetric matrix");
        checkTest(test, testID++);
    }
    testID = 0;
    for (const test of tests.general) {
        checkTest(test, testID++);
    }
})();

describe('Linear solvers (sparse square matrices)', () => {
    describe.each(tests.posDef)('Symmetric positive definite matrices %#', (testCase: TestCase) => {
        expect(testCase.m.isSymmetric()).toBeTruthy();
        describe('Iterative', () => {
            test('ConjGrad', () => {
                let result = ConjugateGradients.solve(testCase.m, testCase.rhs, 20);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(SmallTolerance);
            });
            const initialGuess = Vector.add(testCase.exactSolution, Vector.generate(testCase.exactSolution.size(), (i: number) => { return 0.01 * (Math.random() - 0.5) }));
            test('gauss-zeidel', () => {
                let result = gaussSeidel.solve(testCase.m, testCase.rhs, 120, SmallTolerance, initialGuess);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(Tolerance);
            });
            test('jacobi', () => {
                let result = jacobi.solve(testCase.m, testCase.rhs, 120, SmallTolerance, initialGuess);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(Tolerance);
            });
        });
    });

    describe.each(tests.general)('General matrices %#', (testCase: TestCase) => {
        describe('Iterative', () => {
            const initialGuess = Vector.add(testCase.exactSolution, Vector.generate(testCase.exactSolution.size(), (i: number) => { return (Math.random() - 0.5) }));
            test('gauss-zeidel', () => {
                let result = gaussSeidel.solve(testCase.m, testCase.rhs, 120, SmallTolerance, initialGuess);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(Tolerance);
            });
            test('jacobi', () => {
                let result = jacobi.solve(testCase.m, testCase.rhs, 120, SmallTolerance, initialGuess);
                expect(Vector.lInfDistance(result, testCase.exactSolution)).toBeLessThanOrEqual(Tolerance);
            });
        });
    });
});