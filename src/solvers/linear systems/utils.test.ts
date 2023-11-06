import Matrix from '../../denseMatrix'
import { SmallTolerance } from '../../utils';
import { applyGivensFromLeft, applyGivensFromRight, applyHouseholderFromLeft, applyHouseholderFromRight, applyTransposeGivensFromRight, calcHouseholderVectorCol, calcHouseholderVectorRow, givens, makeGivensMatrix, makeHessenberg } from './eigenvalues';

import { hilbertMatrix, inverseHilbertMatrix } from './utils';

describe('Zeroing methods', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    test('Givens rotations', () => {
        let m = A.clone();
        let m2 = A.clone();
        for (let col = 0; col < m.numCols(); ++col) {
            for (let row = m.numRows() - 1; row > col; --row) {
                let i = row;
                let j = col;
                //const a = m.get(j, 0);
                //const b = m.get(i, 0);
                let givensCoeffs = givens(m.get(j, col), m.get(i, col));
                //console.log(`M[${row}] before: ${m.toString()}`);
                //console.log(`a*c-b*s ${a * givensCoeffs.c - b * givensCoeffs.s}, a*s+b*c${a * givensCoeffs.s + b * givensCoeffs.c}`)
                //console.log(`givensCoeffs a:${a}, b:${b}, c:${givensCoeffs.c}, s:${givensCoeffs.s}, r:${givensCoeffs.r}`);
                applyGivensFromLeft(m, givensCoeffs, i, j);
                //console.log(`M[${row}] after: ${m.toString()}`);
                expect(m.get(j, col)).toBeCloseTo(givensCoeffs.r);
                expect(m.get(i, col)).toBeCloseTo(0);
                m2 = Matrix.mul(makeGivensMatrix(givensCoeffs, A.numRows(), i, j), m2);
                //console.log(`M2[${row}] after: ${m2.toString()}`);
                expect(m2.get(j, col)).toBeCloseTo(givensCoeffs.r);
                expect(m2.get(i, col)).toBeCloseTo(0);
            }
        }
        m = A.clone();
        m2 = A.clone();
        // console.log("cols");
        for (let row = 0; row < m.numRows(); ++row) {
            for (let col = m.numCols() - 1; col > row; --col) {
                let i = col;
                let j = row;
                let givensCoeffs = givens(m.get(row, j), m.get(row, i));
                // todo: remove
                //console.log(`M[${col}] before: ${m.toString()}`);
                //console.log(`givensCoeffs a:${m.get(0, j)}, b:${m.get(0, i)}, c:${givensCoeffs.c}, s:${givensCoeffs.s}, r:${givensCoeffs.r}`);

                applyTransposeGivensFromRight(m, givensCoeffs, i, j);
                //console.log(`M[${col}] after: ${m.toString()}`);
                expect(m.get(row, j)).toBeCloseTo(givensCoeffs.r);
                expect(m.get(row, i)).toBeCloseTo(0);
                m2 = Matrix.mul(m2, makeGivensMatrix(givensCoeffs, A.numRows(), i, j).transposeInPlace());
                //console.log(`M2[${col}] after: ${m2.toString()}`);
                expect(m2.get(row, j)).toBeCloseTo(givensCoeffs.r);
                expect(m2.get(row, i)).toBeCloseTo(0);
            }
        }
    });
    //todo: test symmetric matrix
    test.skip('Householder reflections', () => {
        let m = A.clone();
        let m2 = A.clone();
        for (let col = 0; col + 1 < m.numCols(); ++col) {
            let v = calcHouseholderVectorCol(m, col, col);
            /*applyHouseholderFromLeft(m, v);
            let householderMat = makeHouseholderMatrix(v, col, col);
            m2 = Matrix.mul(householderMat, m2);*/
            for (let row = col + 1; row < m.numRows(); ++row) {
                expect(m.get(row, col)).toBeCloseTo(0);
                expect(m2.get(row, col)).toBeCloseTo(0);
            }
        }
        m = A.clone();
        m2 = A.clone();
        for (let row = 0; row + 1 < m.numRows(); ++row) {
            let v = calcHouseholderVectorRow(m, row, row);
            /*applyHouseholderFromLeft(m, v);
            let householderMat = makeHouseholderMatrix(v, row, row);
            m2 = Matrix.mul(m2, householderMat.transposeInPlace());
            */
            for (let col = row + 1; col < m.numCols(); ++col) {
                expect(m.get(row, col)).toBeCloseTo(0);
                expect(m2.get(row, col)).toBeCloseTo(0);
            }
        }
        throw new Error("Not implemented");

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
});

describe('Hessenberg', () => {
    test('Tridiagonal', () => {
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
        expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, H), Q.transpose()), A)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(H, expectedH)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(Q, expectedQ)).toBeLessThan(SmallTolerance);
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
        console.log(Matrix.mul(Matrix.mul(Q, A), Q.transpose()).toString());
        expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, A), Q.transpose()), H)).toBeLessThan(SmallTolerance);

    });
});

describe('Generators', () => {
    test('Hilberth matrix', () => {
        const hilbert = new Matrix(
            [
                1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
                1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6,
                1 / 3, 1 / 4, 1 / 5, 1 / 6, 1 / 7,
                1 / 4, 1 / 5, 1 / 6, 1 / 7, 1 / 8,
                1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9], 5, 5);
        const invHilbert = new Matrix(
            [
                25, -300, 1050, -1400, 630,
                -300, 4800, -18900, 26880, -12600,
                1050, -18900, 79380, -117600, 56700,
                -1400, 26880, -117600, 179200, -88200,
                630, -12600, 56700, -88200, 44100], 5, 5);
        expect(Matrix.lInfDistance(inverseHilbertMatrix(5), invHilbert)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(hilbertMatrix(5), hilbert)).toBeLessThan(SmallTolerance);
    });
});
