import Matrix from '../../denseMatrix'
import { applyGivensFromLeft, applyGivensFromRight, applyTransposeGivensFromRight, calcHouseholderVectorCol, calcHouseholderVectorRow, givens, makeGivensMatrix, makeHessenberg } from './eigenvalues';

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
        console.log("cols");
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
    });
});

test('Hessenberg', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    let Q: Matrix = Matrix.empty(4, 4);
    let H = makeHessenberg(A, Q);
    expect(H.isHessenberg()).toBeTruthy();
    expect(Q.isOrthogonal()).toBeTruthy();
    expect(Q.isSymmetric()).toBeTruthy();
    expect(Matrix.mul(Q, H));

})