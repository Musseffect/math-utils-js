import Matrix from '../../denseMatrix'
import { applyGivensFromLeft, applyGivensFromRight, applyTransposeGivensFromRight, givens, makeGivensMatrix } from './eigenvalues';

describe('Zeroing methods', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    test('Givens rotations', () => {
        let m = A.clone();
        for (let row = m.numRows() - 1; row > 0; row--) {
            let i = row;
            let j = 0;
            //const a = m.get(j, 0);
            //const b = m.get(i, 0);
            let givensCoeffs = givens(m.get(j, 0), m.get(i, 0));
            //console.log(`M[${row}] before: ${m.toString()}`);
            //console.log(`a*c-b*s ${a * givensCoeffs.c - b * givensCoeffs.s}, a*s+b*c${a * givensCoeffs.s + b * givensCoeffs.c}`)
            //console.log(`givensCoeffs a:${a}, b:${b}, c:${givensCoeffs.c}, s:${givensCoeffs.s}, r:${givensCoeffs.r}`);
            let m2 = m.clone();
            applyGivensFromLeft(m, givensCoeffs, i, j);
            //console.log(`M[${row}] after: ${m.toString()}`);
            expect(m.get(j, 0)).toBeCloseTo(givensCoeffs.r);
            expect(m.get(i, 0)).toBeCloseTo(0);
            m2 = Matrix.mul(makeGivensMatrix(givensCoeffs, A.numRows(), i, j), m2);
            //console.log(`M2[${row}] after: ${m2.toString()}`);
            expect(m2.get(j, 0)).toBeCloseTo(givensCoeffs.r);
            expect(m2.get(i, 0)).toBeCloseTo(0);
        }
        m = A.clone();
        console.log("cols");
        for (let col = m.numCols() - 1; col > 0; col--) {
            let m2 = m.clone();
            let i = col;
            let j = 0;
            let givensCoeffs = givens(m.get(0, j), m.get(0, i));
            // todo: remove
            //console.log(`M[${col}] before: ${m.toString()}`);
            //console.log(`givensCoeffs a:${m.get(0, j)}, b:${m.get(0, i)}, c:${givensCoeffs.c}, s:${givensCoeffs.s}, r:${givensCoeffs.r}`);

            applyTransposeGivensFromRight(m, givensCoeffs, i, j);
            //console.log(`M[${col}] after: ${m.toString()}`);
            expect(m.get(0, j)).toBeCloseTo(givensCoeffs.r);
            expect(m.get(0, i)).toBeCloseTo(0);
            m2 = Matrix.mul(m2, makeGivensMatrix(givensCoeffs, A.numRows(), i, j).transposeInPlace());
            //console.log(`M2[${col}] after: ${m2.toString()}`);
            expect(m2.get(0, j)).toBeCloseTo(givensCoeffs.r);
            expect(m2.get(0, i)).toBeCloseTo(0);
        }
    });
    //todo: test symmetric matrix
    test.skip('Householder reflections', () => {
        throw new Error("Not implemented");
    });
});