
// todo: solvers for rect matrices
// todo: eigenvalue solvers

import Matrix from "../../denseMatrix";
import { Tolerance, SmallTolerance } from "../../utils";
import Vector from "../../vector";
import * as linSolvers from "./exports";


describe.skip('Linear solvers (dense square matrices)', () => {
    interface testData {
        m: Matrix;
        rhs: Vector;
        exactSolution: Vector;
    };
    let testExamples: testData[] = [];

    const mat4 = Matrix.empty(4, 4);
    let rhs = new Vector([6, 25, -11, 15]);
    mat4.set(0, 0, 10);
    mat4.set(0, 1, -1);
    mat4.set(0, 2, 2);

    mat4.set(1, 0, -1);
    mat4.set(1, 1, 11);
    mat4.set(1, 2, -1);
    mat4.set(1, 3, 3);

    mat4.set(2, 0, 2);
    mat4.set(2, 1, -1);
    mat4.set(2, 2, 10);
    mat4.set(2, 3, -1);

    mat4.set(3, 1, 3);
    mat4.set(3, 2, -1);
    mat4.set(3, 3, 8);
    let exactSolution = new Vector([1, 2, -1, 1]);
    testExamples.push({ m: mat4, rhs, exactSolution });
    const choleskyMat = Matrix.empty(3, 3);
    choleskyMat.set(0, 0, 4);
    choleskyMat.set(0, 1, 12);
    choleskyMat.set(0, 2, -16);

    choleskyMat.set(1, 0, 12);
    choleskyMat.set(1, 1, 37);
    choleskyMat.set(1, 2, -43);

    choleskyMat.set(2, 0, -16);
    choleskyMat.set(2, 1, -43);
    choleskyMat.set(2, 2, 98);
    exactSolution = new Vector([1, 2, -3]);
    const choleskyRhs = Matrix.postMulVec(choleskyMat, exactSolution);

    testExamples.push({ m: choleskyMat, rhs: choleskyRhs, exactSolution });

    describe('Symmetric positive definite', ()=> {

        let Q = new Matrix([1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 3, 3, 1, 2, 3, 4], 4, 4);
        let b = new Vector([-4, -7, -9, -10]);
        let x = new Vector([-1, -1, -1, -1]);
        expect(Q.isSymmetric()).toBeTruthy();
        test('ConjGrad', ()=> {
            let solver = new linSolvers.CG();
            linSolvers.CG.solve(Q, b);
        });
        test('LL', ()=> {
            let solver = new linSolvers.LLT();
            expect(solver.factorize(Q)).not.toThrow();
            expect(Vector.sub(solver.solve(b), x).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
        });
        test('LDL', ()=> {
            let solver = new linSolvers.LDLT(SmallTolerance);
            expect(solver.factorize(Q)).not.toThrow();
            expect(Vector.sub(solver.solve(b), x).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
        });
        test('PPLU', ()=>{
            let solve = new linSolvers.PartialPivLU(Q, SmallTolerance);
            //solve.factorize();

        });
        test('FPLU', ()=>{

        });
        test('gauss-zeidel', () => {
            let result = linSolvers.gaussSeidel.solve(Q, b, 30, SmallTolerance);
            expect(Vector.sub(result, x).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
        });
        test('jacobi', () => {
            let result = linSolvers.jacobi.solve(Q, b, 30, SmallTolerance);
            expect(Vector.sub(result, x).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
        });
        test('sor', () => {
            let result = linSolvers.sor.solve(Q, b, 30, 1.0, SmallTolerance);
            expect(Vector.sub(result, x).lInfNorm()).toBeLessThanOrEqual(SmallTolerance);
        });
    });

    describe('Iterative', () => {
        test('gauss-zeidel', () => {
        });
        test('jacobi', () => {

        });
        test('sor', () => {

        });
        test('ConjGrad', () => {
        });
    });

    describe('Decompositions', () => {
        test('PartialPivLU', () => {
            for (let testExample of testExamples) {
                expect(Vector.near(linSolvers.PartialPivLU.solve(testExample.m.clone(), testExample.rhs.clone(), SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();
                // test decomposition
                let luSolver = new linSolvers.PartialPivLU(testExample.m);
                let PT = luSolver.P().transpose();
                let PTLU = Matrix.mul(PT, Matrix.mul(luSolver.L().toMatrix(), luSolver.U().toMatrix()));
                expect(Math.abs(PT.determinantNaive())).toBeCloseTo(1);
                expect(Matrix.near(PTLU, testExample.m));

                expect(Matrix.near(linSolvers.PartialPivLU.solveMatrix(testExample.m.clone(), Matrix.identity(testExample.m.height())), testExample.m.inverseNaive()));
            }
        });

        test('FullPivLU', () => {
            for (let testExample of testExamples) {
                expect(Vector.near(linSolvers.FullPivLU.solve(testExample.m.clone(), testExample.rhs.clone(), SmallTolerance), testExample.exactSolution, Tolerance)).toBeTruthy();
                // test decomposition
                let luSolver = new linSolvers.FullPivLU(testExample.m);
                let PT = luSolver.P().transpose();
                let QT = luSolver.Q().transpose();
                let PTLUQT = Matrix.mul(Matrix.mul(PT, Matrix.mul(luSolver.L().toMatrix(), luSolver.U().toMatrix())), QT);
                expect(Math.abs(PT.determinantNaive())).toBeCloseTo(1);
                expect(Math.abs(QT.determinantNaive())).toBeCloseTo(1);
                expect(Matrix.near(PTLUQT, testExample.m));

                expect(Matrix.near(linSolvers.PartialPivLU.solveMatrix(testExample.m.clone(), Matrix.identity(testExample.m.height())), testExample.m.inverseNaive()));
            }
        });
        test('LL', () => {

        });
        test('LDL', () => {

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