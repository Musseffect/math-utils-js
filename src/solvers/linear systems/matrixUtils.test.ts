import Matrix from '../../denseMatrix'
import { SmallTolerance, StopWatch, sign } from '../../utils';
import { givens, applyGivensFromLeft, applyTransposeGivensFromRight, makeGivensMatrix, applyGivensFromRight, applyTransposeGivensFromLeft } from './givensRotation';
import { applyHouseholderFromLeft, applyHouseholderFromRight, calcHouseholderVectorCol, calcHouseholderVectorRow, makeHouseholder } from './hausholderReflection';
import { makeHessenberg, makeTridiagonal, makeTridiagonalAlt } from './hessenbergMatrix';

import { MatrixGenerator } from '../../dense/matrixGenerator';

import fs from 'fs';
//import * as v8Profiler from 'v8-profiler-next';
import JSGenerator from '../../random/js';

// QR decomposition
describe.skip('Upper triangular zeroing', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    test('Givens rotations: implicit', () => {
        let R = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let col = 0; col < R.numCols(); ++col) {
            for (let row = R.numRows() - 1; row > col; --row) {
                let i = row;
                let j = col;
                let givensCoeffs = givens(R.get(j, col), R.get(i, col));
                applyGivensFromLeft(R, givensCoeffs, i, j);
                applyTransposeGivensFromRight(Q, givensCoeffs, i, j);
                expect(R.get(j, col)).toBeCloseTo(givensCoeffs.r);
                expect(R.get(i, col)).toBeCloseTo(0);
            }
        }
        expect(R.isTriangular(true)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(Q, R))).toBeLessThan(SmallTolerance);
    });
    // QR with givens rotations
    test('Givens rotation: explicit', () => {
        let R = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let col = 0; col < R.numCols(); ++col) {
            for (let row = R.numRows() - 1; row > col; --row) {
                let i = row;
                let j = col;
                let givensCoeffs = givens(R.get(j, col), R.get(i, col));
                let Q_k = makeGivensMatrix(givensCoeffs, R.numRows(), i, j);
                expect(Q_k.isOrthogonal()).toBeTruthy();
                Q = Matrix.mul(Q, Q_k.transpose());
                R = Matrix.mul(Q_k, R);
                expect(R.get(j, col)).toBeCloseTo(givensCoeffs.r);
                expect(R.get(i, col)).toBeCloseTo(0);
            }
        }
        expect(R.isTriangular(true)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(Q, R))).toBeLessThan(SmallTolerance);
    })
    //todo: test symmetric matrix
    test('Householder reflections: implicit', () => {
        let R = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let col = 0; col + 1 < R.numCols(); ++col) {
            let xNorm = 0.0;
            for (let row = col; row < R.numRows(); ++row)
                xNorm += Math.pow(R.get(row, col), 2);
            xNorm = -sign(R.get(col, col)) * Math.sqrt(xNorm);
            let v = calcHouseholderVectorCol(R, col, col);
            applyHouseholderFromLeft(v, R, col);
            applyHouseholderFromRight(v, Q, col);
            expect(R.get(col, col)).toBeCloseTo(xNorm);
            for (let row = col + 1; row < R.numRows(); ++row)
                expect(R.get(row, col)).toBeCloseTo(0);
        }
        expect(R.isTriangular(true)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(Q, R))).toBeLessThan(SmallTolerance);
        /*m = A.clone();
        m2 = A.clone();
        for (let row = 0; row + 1 < m.numRows(); ++row) {
            let v = calcHouseholderVectorRow(m, row, row);
            applyHouseholderFromLeft(m, v);
            let householderMat = makeHouseholderMatrix(v, row, row);
            m2 = Matrix.mul(m2, householderMat.transposeInPlace());
            for (let col = row + 1; col < m.numCols(); ++col) {
                expect(m.get(row, col)).toBeCloseTo(0);
                expect(m2.get(row, col)).toBeCloseTo(0);
            }
        }
        throw new Error("Not implemented");
            */

        // test symmetric householder by generating hessenberg matrix QAQT = H
        /*let A = new Matrix([
            4, 1, -2, 2,
            1, 2, 0, 1,
            -2, 0, 3, -2,
            2, 1, -2, -1], 4, 4);
        let Q = Matrix.identity(A.numRows());
        let H = A.clone();
        for (let iter = 0; iter + 2 < A.numCols(); ++iter) {
            let v = calcHouseholderVectorRow(m, iter + 2, iter);
            applyHessenbergFromLeft(v, Q, iter + 2);
            applyHessenbergFromLeft(v, H, iter + 2);
            applyHessenbergFromRight(v, H, iter + 2);
            // check matrices
        }
        expect(H.isHessenberg()).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(Matrix.mul(Q, Matrix.mul(A, Q.transpose())), H)).toBeLessThan(SmallTolerance);
    */
    });
    test('Householder rotation: explicit', () => {
        let R = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let col = 0; col + 1 < R.numCols(); ++col) {
            let xNorm = 0.0;
            for (let row = col; row < R.numRows(); ++row)
                xNorm += Math.pow(R.get(row, col), 2);
            xNorm = -sign(R.get(col, col)) * Math.sqrt(xNorm);
            let v = calcHouseholderVectorCol(R, col, col);
            let Q_k = makeHouseholder(v, R.numRows());
            expect(Q_k.isOrthogonal()).toBeTruthy();
            Q = Matrix.mul(Q, Q_k);
            R = Matrix.mul(Q_k, R);
            expect(R.get(col, col)).toBeCloseTo(xNorm);
            for (let row = col + 1; row < R.numRows(); ++row)
                expect(R.get(row, col)).toBeCloseTo(0);
        }
        expect(R.isTriangular(true)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(Q, R))).toBeLessThan(SmallTolerance);
    });
});

// LQ decomposition
describe.skip('Lower triangular zeroing', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    test.skip('Givens rotation: implicit', () => {
        let L = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let row = 0; row < L.numRows(); ++row) {
            for (let col = L.numCols() - 1; col > row; --col) {
                let i = col;
                let j = row;
                let givensCoeffs = givens(L.get(row, j), L.get(row, i));
                applyGivensFromRight(L, givensCoeffs, i, j);
                applyTransposeGivensFromLeft(Q, givensCoeffs, i, j);
                expect(L.get(row, j)).toBeCloseTo(givensCoeffs.r);
                expect(L.get(row, i)).toBeCloseTo(0);
            }
        }
        expect(L.isTriangular(false)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(L, Q))).toBeLessThan(SmallTolerance);
    });
    test.skip('Givens rotation: explicit', () => {
        let L = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let row = 0; row < L.numRows(); ++row) {
            throw new Error("Not implemented");
        }
        expect(L.isTriangular(false)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(L, Q))).toBeLessThan(SmallTolerance);
    });
    test.skip('Householder rotation: implicit', () => {
        let L = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let row = 0; row < L.numRows(); ++row) {

            throw new Error("Not implemented");
        }
        expect(L.isTriangular(false)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(L, Q))).toBeLessThan(SmallTolerance);
    });
    test.skip('Householder rotation: explicit', () => {
        let L = A.clone();
        let Q = Matrix.identity(A.numRows());
        for (let row = 0; row < L.numRows(); ++row) {
            throw new Error("Not implemented");
        }
        expect(L.isTriangular(false)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(A, Matrix.mul(L, Q))).toBeLessThan(SmallTolerance);
    });
});

describe.skip('Upper hessenberg zeroing', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    // QAQT = H;
    test('Givens rotations: explicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let col = 0; col + 2 < H.numCols(); ++col) {
            for (let row = H.numRows() - 1; row > col + 1; --row) {
                const i = row;
                const j = col + 1;
                let givensCoeffs = givens(H.get(j, col), H.get(i, col));
                let Q_k = makeGivensMatrix(givensCoeffs, H.numRows(), i, j);
                H = Matrix.mul(Matrix.mul(Q_k, H), Q_k.transpose());
                Q = Matrix.mul(Q_k, Q);
                expect(H.get(j, col)).toBeCloseTo(givensCoeffs.r);
                expect(H.get(i, col)).toBeCloseTo(0);
            }
        }
        expect(H.isHessenberg(true));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
    test('Givens rotations: implicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let col = 0; col + 2 < H.numCols(); ++col) {
            for (let row = H.numRows() - 1; row > col + 1; --row) {
                const i = row;
                const j = col + 1;
                let givensCoeffs = givens(H.get(j, col), H.get(i, col));
                applyGivensFromLeft(H, givensCoeffs, i, j);
                applyTransposeGivensFromRight(H, givensCoeffs, i, j);
                applyGivensFromLeft(Q, givensCoeffs, i, j);
                expect(H.get(j, col)).toBeCloseTo(givensCoeffs.r);
                expect(H.get(i, col)).toBeCloseTo(0);
            }
        }
        expect(H.isHessenberg(true));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
    test('Householder reflections:explicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let col = 0; col + 2 < H.numCols(); ++col) {
            let xNorm = 0.0;
            for (let row = col + 1; row < H.numRows(); ++row)
                xNorm += Math.pow(H.get(row, col), 2);
            xNorm = -sign(H.get(col + 1, col)) * Math.sqrt(xNorm);
            let v = calcHouseholderVectorCol(H, col + 1, col);
            let Q_k = makeHouseholder(v, H.numRows());
            expect(Q_k.isOrthogonal()).toBeTruthy();
            expect(Q_k.isSymmetric()).toBeTruthy();
            Q = Matrix.mul(Q_k, Q);
            H = Matrix.mul(Matrix.mul(Q_k, H), Q_k);
            expect(H.get(col + 1, col)).toBeCloseTo(xNorm);
            for (let row = col + 2; row < H.numRows(); ++row)
                expect(H.get(row, col)).toBeCloseTo(0);
        }
        expect(H.isHessenberg(true));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
    test('Householder reflections:implicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let col = 0; col + 2 < H.numCols(); ++col) {
            let xNorm = 0.0;
            for (let row = col + 1; row < H.numRows(); ++row)
                xNorm += Math.pow(H.get(row, col), 2);
            xNorm = -sign(H.get(col + 1, col)) * Math.sqrt(xNorm);
            let v = calcHouseholderVectorCol(H, col + 1, col);
            applyHouseholderFromLeft(v, H, col + 1);
            applyHouseholderFromRight(v, H, col + 1);
            applyHouseholderFromLeft(v, Q, col + 1);
            expect(H.get(col + 1, col)).toBeCloseTo(xNorm);
            for (let row = col + 2; row < H.numRows(); ++row)
                expect(H.get(row, col)).toBeCloseTo(0);
        }
        expect(H.isHessenberg(true));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
});

describe.skip('Lower hessenberg zeroing', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    test('Givens rotations: explicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let row = 0; row + 2 < H.numRows(); ++row) {
            for (let col = H.numCols() - 1; col > row + 1; --col) {
                const i = col;
                const j = row + 1;
                let givensCoeffs = givens(H.get(row, j), H.get(row, i));
                let Q_k = makeGivensMatrix(givensCoeffs, H.numRows(), i, j);
                H = Matrix.mul(Matrix.mul(Q_k, H), Q_k.transpose());
                Q = Matrix.mul(Q_k, Q);
                expect(H.get(row, j)).toBeCloseTo(givensCoeffs.r);
                expect(H.get(row, i)).toBeCloseTo(0);
            }
        }
        expect(H.isHessenberg(false));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
    test('Givens rotations: implicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let row = 0; row + 2 < H.numRows(); ++row) {
            for (let col = H.numCols() - 1; col > row + 1; --col) {
                const i = col;
                const j = row + 1;
                let givensCoeffs = givens(H.get(row, j), H.get(row, i));
                applyGivensFromLeft(H, givensCoeffs, i, j);
                applyTransposeGivensFromRight(H, givensCoeffs, i, j);
                applyGivensFromLeft(Q, givensCoeffs, i, j);
                expect(H.get(row, j)).toBeCloseTo(givensCoeffs.r);
                expect(H.get(row, i)).toBeCloseTo(0);
            }
        }
        expect(H.isHessenberg(false));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
    test('Householder reflections:explicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let row = 0; row + 2 < H.numRows(); ++row) {
            let xNorm = 0.0;
            for (let col = row + 1; col < H.numRows(); ++col)
                xNorm += Math.pow(H.get(row, col), 2);
            xNorm = -sign(H.get(row, row + 1)) * Math.sqrt(xNorm);
            let v = calcHouseholderVectorRow(H, row, row + 1);
            let Q_k = makeHouseholder(v, H.numCols());
            expect(Q_k.isOrthogonal()).toBeTruthy();
            expect(Q_k.isSymmetric()).toBeTruthy();
            Q = Matrix.mul(Q_k, Q);
            H = Matrix.mul(Matrix.mul(Q_k, H), Q_k);
            expect(H.get(row, row + 1)).toBeCloseTo(xNorm);
            for (let col = row + 2; col < H.numRows(); ++col)
                expect(H.get(row, col)).toBeCloseTo(0);
        }
        expect(H.isHessenberg(false));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
    test('Householder reflections:implicit', () => {
        let H = A.clone();
        let Q = Matrix.identity(H.numRows());
        for (let row = 0; row + 2 < H.numRows(); ++row) {
            let xNorm = 0.0;
            for (let col = row + 1; col < H.numRows(); ++col)
                xNorm += Math.pow(H.get(row, col), 2);
            xNorm = -sign(H.get(row, row + 1)) * Math.sqrt(xNorm);
            let v = calcHouseholderVectorRow(H, row, row + 1);
            applyHouseholderFromLeft(v, H, row + 1);
            applyHouseholderFromRight(v, H, row + 1);
            applyHouseholderFromLeft(v, Q, row + 1);
            expect(H.get(row, row + 1)).toBeCloseTo(xNorm);
            for (let col = row + 2; col < H.numRows(); ++col)
                expect(H.get(row, col)).toBeCloseTo(0);
        }
        expect(H.isHessenberg(false));
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(H, Matrix.mul(Matrix.mul(Q, A), Q.transpose()))).toBeLessThan(SmallTolerance);
    });
});

//v8Profiler.setGenerateType(1);
const title = 'Hessenberg-performance';

function startPerformanceProfiling() {
    /*v8Profiler.startProfiling(title, true);
    afterAll(() => {
        const profile = v8Profiler.stopProfiling(title);
        profile.export(function (error, result: any) {
            // if it doesn't have the extension .cpuprofile then
            // chrome's profiler tool won't like it.
            // examine the profile:
            //   Navigate to chrome://inspect
            //   Click Open dedicated DevTools for Node
            //   Select the profiler tab
            //   Load your file
            fs.writeFileSync(`${title}.cpuprofile`, result);
            profile.delete();
        });
    });*/
}

describe.skip('Hessenberg performance', () => {
    const numRepetitions = 20;
    let matrices: Matrix[] = [];
    for (const size of [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        matrices.push(Matrix.random(size, size));

    // startPerformanceProfiling();

    test.skip('Hessenberg partial default', () => {
        let stopWatch = new StopWatch();
        let time: number[] = [];
        for (const A of matrices) {
            let H: Matrix = A;
            stopWatch.reset();
            for (let i = 0; i < numRepetitions; ++i)
                H = makeHessenberg(A, undefined);
            time.push(stopWatch.elapsed() / numRepetitions);
            expect(H.isHessenberg(true)).toBeTruthy();
        }
        console.log(`Hessenberg partial default: ${time}`);
    });
    test('Full Hessenberg default', () => {
        let stopWatch = new StopWatch();
        let time: number[] = [];
        for (const A of matrices) {
            let Q: Matrix = Matrix.empty(A.width(), A.width());
            let H: Matrix = A;
            stopWatch.reset();
            for (let i = 0; i < numRepetitions; ++i)
                H = makeHessenberg(A, Q);
            time.push(stopWatch.elapsed() / numRepetitions);
            expect(H.isHessenberg(true)).toBeTruthy();
            expect(Q.isOrthogonal()).toBeTruthy();
            expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, A), Q.transpose()), H)).toBeLessThan(SmallTolerance);
        }
        console.log(`Full hessenberg default: ${time}`);
    });
});

test.skip('Triangular alt', () => {
    let S: Matrix = new Matrix([
        4, 1, -2, 2,
        1, 2, 0, 1,
        -2, 0, 3, -2,
        2, 1, -2, -1], 4, 4);
    let expectedQ = new Matrix([
        1, 0, 0, 0,
        0, -1 / 3, 2 / 3, -2 / 3,
        0, 2 / 15, -2 / 3, -11 / 15,
        0, -14 / 15, -1 / 3, 2 / 15
    ], 4, 4);
    let expectedT = new Matrix([
        4, -3, 0, 0,
        -3, 10 / 3, -5 / 3, 0,
        0, -5 / 3, -33 / 25, 68 / 75,
        0, 0, 68 / 75, 149 / 75], 4, 4);
    let Q: Matrix = Matrix.empty(4, 4);
    let T1 = makeTridiagonalAlt(S);
    let T = makeTridiagonalAlt(S, Q);
    let TOld = makeTridiagonal(S);
    console.log(`Expected: ${expectedT.toString()}`);
    console.log(`Told: ${TOld.toString()}`);
    console.log(`TAlt1: ${T1.toString()}`);
    console.log(`TAlt2: ${T.toString()}`);
    expect(Matrix.lInfDistance(T, T1)).toBeLessThan(SmallTolerance);
    expect(Matrix.lInfDistance(TOld, T)).toBeLessThan(SmallTolerance);
    expect(T.isHessenberg(true)).toBeTruthy();
    expect(T.isTridiagonal()).toBeTruthy();
    expect(Q.isOrthogonal()).toBeTruthy();
    expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, S), Q.transpose()), T)).toBeLessThan(SmallTolerance);
    expect(Matrix.lInfDistance(T, expectedT)).toBeLessThan(SmallTolerance);
    expect(Matrix.lInfDistance(Q, expectedQ)).toBeLessThan(SmallTolerance);
});

describe('Triangular performance', () => {
    const generator = new MatrixGenerator(new JSGenerator());
    const numRepetitions = 10;
    let matrices: Matrix[] = [];
    for (const size of [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        matrices.push(generator.randomSymmetric(size));

    startPerformanceProfiling();

    test('Partial Triangular alt', () => {
        //let stopWatch = new StopWatch();
        //let time: number[] = [];
        for (const S of matrices) {
            //stopWatch.reset();
            let T: Matrix = S;
            for (let i = 0; i < numRepetitions; ++i)
                T = makeTridiagonalAlt(S, undefined);
            //time.push(stopWatch.elapsed());
            //expect(T.isTridiagonal()).toBeTruthy();
        }
        //console.log(`Partial Triangular alt' Q: ${time}`);
    });
    test('Full Triangular alt', () => {
        //let stopWatch = new StopWatch();
        //let time: number[] = [];
        for (const S of matrices) {
            let Q: Matrix = Matrix.empty(S.width(), S.width());
            //stopWatch.reset();
            let T: Matrix = S;
            for (let i = 0; i < numRepetitions; ++i) {
                T = makeTridiagonalAlt(S, Q);
            }
            //time.push(stopWatch.elapsed());
            //expect(T.isTridiagonal()).toBeTruthy();
            //expect(Q.isOrthogonal()).toBeTruthy();
            //expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, S), Q.transpose()), T)).toBeLessThan(SmallTolerance);
        }
        //console.log(`Full Triangular alt Q: ${time}`);
    });
    test('Partial Triangular default', () => {
        //let stopWatch = new StopWatch();
        //let time: number[] = [];
        for (const S of matrices) {
            //stopWatch.reset();
            let T: Matrix = S;
            for (let i = 0; i < numRepetitions; ++i)
                T = makeTridiagonal(S, undefined);
            //time.push(stopWatch.elapsed());
            //expect(T.isTridiagonal()).toBeTruthy();
        }
        //console.log(`Partial Triangular default: ${time}`);
    });
    test('Full Triangular default', () => {
        //let stopWatch = new StopWatch();
        //let time: number[] = [];
        for (const S of matrices) {
            let Q: Matrix = Matrix.empty(S.width(), S.width());
            //stopWatch.reset();
            let T: Matrix = S;
            for (let i = 0; i < numRepetitions; ++i)
                T = makeTridiagonal(S, Q);
            //time.push(stopWatch.elapsed());
            //expect(T.isTridiagonal()).toBeTruthy();
            //expect(Q.isOrthogonal()).toBeTruthy();
            //expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, S), Q.transpose()), T)).toBeLessThan(SmallTolerance);

        }
        //console.log(`Full Triangular default: ${time}`);
    });
});

describe.skip('Hessenberg form decomposition', () => {
    test.skip('Tridiagonal', () => {
        let A: Matrix = new Matrix([
            4, 1, -2, 2,
            1, 2, 0, 1,
            -2, 0, 3, -2,
            2, 1, -2, -1], 4, 4);
        let expectedQ = new Matrix([
            1, 0, 0, 0,
            0, -1 / 3, 2 / 3, -2 / 3,
            0, 2 / 15, -2 / 3, -11 / 15,
            0, -14 / 15, -1 / 3, 2 / 15
        ], 4, 4);
        let expectedH = new Matrix([
            4, -3, 0, 0,
            -3, 10 / 3, -5 / 3, 0,
            0, -5 / 3, -33 / 25, 68 / 75,
            0, 0, 68 / 75, 149 / 75], 4, 4);
        let Q: Matrix = Matrix.empty(4, 4);
        let H = makeHessenberg(A, Q);
        expect(H.isHessenberg(true)).toBeTruthy();
        expect(H.isTridiagonal()).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, A), Q.transpose()), H)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(H, expectedH)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(Q, expectedQ)).toBeLessThan(SmallTolerance);
        Q = Matrix.empty(4, 4);
        let T = makeTridiagonal(A, Q);
        expect(T.isTridiagonal()).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(T, expectedH)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(Q, expectedQ)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, A), Q.transpose()), T)).toBeLessThan(SmallTolerance);
    });
    test('Hessenberg', () => {
        let A: Matrix = new Matrix([
            1, 2, 3, 4,
            5, 6, 7, 8,
            3, 4, 2, -2,
            3, 5, 1, 2], 4, 4);
        let Q: Matrix = Matrix.empty(4, 4);
        let H = makeHessenberg(A, Q);
        expect(H.isHessenberg(true)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, A), Q.transpose()), H)).toBeLessThan(SmallTolerance);
    });
});