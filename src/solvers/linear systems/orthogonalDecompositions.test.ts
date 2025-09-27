import Matrix from '../../dense/denseMatrix'
import { SmallTolerance, SmallestTolerance, StopWatch, Tolerance, assert, sign } from '../../utils';
import { givens, applyGivensFromLeft, applyTransposeGivensFromRight, makeGivensMatrix, applyGivensFromRight, applyTransposeGivensFromLeft } from './givensRotation';
import { applyHouseholderFromLeft, applyHouseholderFromRight, calcHouseholderVectorCol, calcHouseholderVectorInplace, calcHouseholderVectorRow, makeHouseholderMatrix } from './hausholderReflection';
import { makeHessenberg, makeTridiagonal, makeTridiagonalAlt } from './hessenbergMatrix';

import { MatrixGenerator } from '../../dense/matrixGenerator';

import fs from 'fs';
//import * as v8Profiler from 'v8-profiler-next';
import JSGenerator from '../../random/js';
import { jacobiRotation } from './jacobiRotation';
import { OrthogonalDecomposition, OrthogonalDecompositionType, OrthogonalDecompositionParams, ZeroingMethod, QRTest } from './qr';
import Vector from '../../dense/vector';
import { PermutationMatrix, PermutationType } from '../../permutationMatrix';

describe('Transformations', () => {
    test('Jacobi', () => {
        let a1 = 2;
        let a2 = 1;
        let b = 11;
        const { c, s } = jacobiRotation(a1, a2, b);
        expect(c * c + s * s).toBeCloseTo(1);
        expect(b * (c * c - s * s) + (a1 - a2) * c * s).toBeCloseTo(0);
        let J = new Matrix([c, s, -s, c], 2, 2);
        let A = new Matrix([a1, b, b, a2], 2, 2);
        let result = Matrix.mul(Matrix.mul(J.transpose(), A), J);
        expect(result.isDiagonal()).toBeTruthy();
    });
    test('Householder', () => {
        const Size = 4;
        let emptyVec = Vector.empty(Size);
        expect(calcHouseholderVectorInplace(emptyVec, 0).l2Norm()).toBeCloseTo(0);
        let nonEmptyVec = new Vector([1, 2, 3, 4]);
        let result = calcHouseholderVectorInplace(nonEmptyVec.getSubVector(1, Size - 1), 0);
        let M = makeHouseholderMatrix(result, 1, Size);
        let reflected = Matrix.postMulVec(M, nonEmptyVec);
        expect(reflected.get(0)).not.toBeCloseTo(0);
        expect(reflected.get(1)).not.toBeCloseTo(0);
        expect(reflected.get(2)).toBeCloseTo(0);
        expect(reflected.get(3)).toBeCloseTo(0);
    });
    test('Givens', () => {
        let a = 0;
        let b = 0;
        const { c, s, r } = givens(a, b);
        expect(r).toBeCloseTo(0);
        expect(c).toBeCloseTo(1);
        expect(s).toBeCloseTo(0);
    })
});

describe('Upper triangular zeroing', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    describe("QR", () => {
        const checkResult = (Q: Matrix, R: Matrix) => {
            expect(R.isTriangular(true)).toBeTruthy();
            expect(Q.isOrthogonal()).toBeTruthy();
            expect(Matrix.lInfDistance(A, Matrix.mul(Q, R))).toBeLessThan(SmallTolerance);
        };
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
            checkResult(Q, R);
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
            checkResult(Q, R);
        })
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
            checkResult(Q, R);
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
                let Q_k = makeHouseholderMatrix(v, R.numRows() - v.size(), R.numRows());
                expect(Q_k.isOrthogonal()).toBeTruthy();
                Q = Matrix.mul(Q, Q_k);
                R = Matrix.mul(Q_k, R);
                expect(R.get(col, col)).toBeCloseTo(xNorm);
                for (let row = col + 1; row < R.numRows(); ++row)
                    expect(R.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, R);
        });
    });
    describe("RQ", () => {
        const checkResult = (Q: Matrix, R: Matrix) => {
            expect(R.isTriangular(true)).toBeTruthy();
            expect(Q.isOrthogonal()).toBeTruthy();
            expect(Matrix.lInfDistance(A, Matrix.mul(R, Q))).toBeLessThan(SmallTolerance);
        };
        test('Givens rotations: implicit', () => {
            let R = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = R.numRows() - 1; row > 0; --row) {
                for (let col = 0; col < row; ++col) {
                    let i = row;
                    let j = col;
                    let givensCoeffs = givens(R.get(row, i), R.get(row, j));
                    applyGivensFromRight(R, givensCoeffs, i, j);
                    applyTransposeGivensFromLeft(Q, givensCoeffs, i, j);
                    expect(R.get(row, i)).toBeCloseTo(givensCoeffs.r);
                    expect(R.get(row, j)).toBeCloseTo(0);
                }
            }
            checkResult(Q, R);
        });
        test('Givens rotations: explicit', () => {
            let R = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = R.numRows() - 1; row > 0; --row) {
                for (let col = 0; col < row; ++col) {
                    let i = row;
                    let j = col;
                    let givensCoeffs = givens(R.get(row, i), R.get(row, j));
                    let Q_k = makeGivensMatrix(givensCoeffs, R.numRows(), i, j);
                    expect(Q_k.isOrthogonal()).toBeTruthy();
                    Q = Matrix.mul(Q_k.transpose(), Q);
                    R = Matrix.mul(R, Q_k);
                    expect(R.get(row, i)).toBeCloseTo(givensCoeffs.r);
                    expect(R.get(row, j)).toBeCloseTo(0);
                }
            }
            checkResult(Q, R);
        });
        test('Householder rotation: implicit', () => {
            let R = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = A.numRows() - 1; row > 0; --row) {
                let xNorm = 0.0;
                for (let col = 0; col <= row; ++col)
                    xNorm += Math.pow(R.get(row, col), 2);
                xNorm = -sign(R.get(row, row)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorRow(R, row, 0, row + 1, row);
                applyHouseholderFromRight(v, R, 0);
                applyHouseholderFromLeft(v, Q, 0);
                expect(R.get(row, row)).toBeCloseTo(xNorm);
                for (let col = 0; col < row; ++col)
                    expect(R.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, R);
        });
        test('Householder rotation: explicit', () => {
            let R = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = A.numRows() - 1; row > 0; --row) {
                let xNorm = 0.0;
                for (let col = 0; col <= row; ++col)
                    xNorm += Math.pow(R.get(row, col), 2);
                xNorm = -sign(R.get(row, row)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorRow(R, row, 0, row + 1, row);
                let Q_k = makeHouseholderMatrix(v, 0, R.numRows());
                expect(Q_k.isOrthogonal()).toBeTruthy();
                Q = Matrix.mul(Q_k, Q);
                R = Matrix.mul(R, Q_k);
                expect(R.get(row, row)).toBeCloseTo(xNorm);
                for (let col = 0; col < row; ++col)
                    expect(R.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, R);
        });
    });
});

describe('Lower triangular zeroing', () => {
    let A: Matrix = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4);
    describe("LQ", () => {
        const checkResult = (Q: Matrix, L: Matrix) => {
            expect(L.isTriangular(false)).toBeTruthy();
            expect(Q.isOrthogonal()).toBeTruthy();
            expect(Matrix.lInfDistance(A, Matrix.mul(L, Q))).toBeLessThan(SmallTolerance);
        };
        test('Givens rotation: implicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = 0; row < L.numRows(); ++row) {
                for (let col = L.numCols() - 1; col > row; --col) {
                    let i = col;
                    let j = row;
                    let givensCoeffs = givens(L.get(row, j), L.get(row, i));
                    applyTransposeGivensFromRight(L, givensCoeffs, i, j);
                    applyGivensFromLeft(Q, givensCoeffs, i, j);
                    expect(L.get(row, j)).toBeCloseTo(givensCoeffs.r);
                    expect(L.get(row, i)).toBeCloseTo(0);
                }
            }
            checkResult(Q, L);
        });
        test('Givens rotation: explicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = 0; row < L.numRows(); ++row) {
                for (let col = L.numCols() - 1; col > row; --col) {
                    let i = col;
                    let j = row;
                    let givensCoeffs = givens(L.get(row, j), L.get(row, i));
                    let Q_k = makeGivensMatrix(givensCoeffs, L.numCols(), i, j);
                    expect(Q_k.isOrthogonal()).toBeTruthy();
                    L = Matrix.mul(L, Q_k.transpose());
                    Q = Matrix.mul(Q_k, Q);
                    expect(L.get(row, j)).toBeCloseTo(givensCoeffs.r);
                    expect(L.get(row, i)).toBeCloseTo(0);
                }
            }
            checkResult(Q, L);
        });
        test('Householder rotation: implicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = 0; row < L.numRows(); ++row) {
                let xNorm = 0.0;
                for (let col = row; col < L.numCols(); ++col)
                    xNorm += Math.pow(L.get(row, col), 2);
                xNorm = -sign(L.get(row, row)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorRow(L, row, row);
                applyHouseholderFromRight(v, L, row);
                applyHouseholderFromLeft(v, Q, row);
                expect(L.get(row, row)).toBeCloseTo(xNorm);
                for (let col = row + 1; col < L.numCols(); ++col)
                    expect(L.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, L);
        });
        test('Householder rotation: explicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let row = 0; row < L.numRows(); ++row) {
                let xNorm = 0.0;
                for (let col = row; col < L.numCols(); ++col)
                    xNorm += Math.pow(L.get(row, col), 2);
                xNorm = -sign(L.get(row, row)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorRow(L, row, row);
                let Q_k = makeHouseholderMatrix(v, L.numRows() - v.size(), L.numRows());
                expect(Q_k.isOrthogonal()).toBeTruthy();
                Q = Matrix.mul(Q_k, Q);
                L = Matrix.mul(L, Q_k);

                expect(L.get(row, row)).toBeCloseTo(xNorm);
                for (let col = row + 1; col < L.numCols(); ++col)
                    expect(L.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, L);
        });
    });
    describe("QL", () => {
        const checkResult = (Q: Matrix, L: Matrix) => {
            expect(L.isTriangular(false)).toBeTruthy();
            expect(Q.isOrthogonal()).toBeTruthy();
            expect(Matrix.lInfDistance(A, Matrix.mul(Q, L))).toBeLessThan(SmallTolerance);
        };
        test('Givens rotations: implicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let col = L.numCols() - 1; col > 0; --col) {
                for (let row = 0; row < col; ++row) {
                    let i = col;
                    let j = row;
                    let givensCoeffs = givens(L.get(i, col), L.get(j, col));
                    applyTransposeGivensFromLeft(L, givensCoeffs, i, j);
                    applyGivensFromRight(Q, givensCoeffs, i, j);
                    expect(L.get(i, col)).toBeCloseTo(givensCoeffs.r);
                    expect(L.get(j, col)).toBeCloseTo(0);
                }
            }
            checkResult(Q, L);
        });
        test('Givens rotations: explicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let col = L.numCols() - 1; col > 0; --col) {
                for (let row = 0; row < col; ++row) {
                    let i = col;
                    let j = row;
                    let givensCoeffs = givens(L.get(i, col), L.get(j, col));
                    let Q_k = makeGivensMatrix(givensCoeffs, L.numCols(), i, j);
                    expect(Q_k.isOrthogonal()).toBeTruthy();
                    L = Matrix.mul(Q_k.transpose(), L);
                    Q = Matrix.mul(Q, Q_k);
                    expect(L.get(i, col)).toBeCloseTo(givensCoeffs.r);
                    expect(L.get(j, col)).toBeCloseTo(0);
                }
            }
            checkResult(Q, L);
        });
        test('Householder rotation: implicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let col = L.numCols() - 1; col > 0; --col) {
                let xNorm = 0.0;
                for (let row = 0; row <= col; ++row)
                    xNorm += Math.pow(L.get(row, col), 2);
                xNorm = -sign(L.get(col, col)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorCol(L, 0, col, col + 1, col);
                applyHouseholderFromLeft(v, L, 0);
                applyHouseholderFromRight(v, Q, 0);

                expect(L.get(col, col)).toBeCloseTo(xNorm);
                for (let row = 0; row < col; ++row)
                    expect(L.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, L);
        });
        test('Householder rotation: explicit', () => {
            let L = A.clone();
            let Q = Matrix.identity(A.numRows());
            for (let col = L.numCols() - 1; col > 0; --col) {
                let xNorm = 0.0;
                for (let row = 0; row <= col; ++row)
                    xNorm += Math.pow(L.get(row, col), 2);
                xNorm = -sign(L.get(col, col)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorCol(L, 0, col, col + 1, col);
                let Q_k = makeHouseholderMatrix(v, 0, L.numRows());
                expect(Q_k.isOrthogonal()).toBeTruthy();
                Q = Matrix.mul(Q, Q_k);
                L = Matrix.mul(Q_k, L);

                expect(L.get(col, col)).toBeCloseTo(xNorm);
                for (let row = 0; row < col; ++row)
                    expect(L.get(row, col)).toBeCloseTo(0);
            }
            checkResult(Q, L);
        });
    });
});

describe('Upper hessenberg zeroing', () => {
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
            let v = calcHouseholderVectorCol(H, col + 1, col, H.numRows() - col - 1);
            let Q_k = makeHouseholderMatrix(v, H.numRows() - v.size(), H.numRows());
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
            let v = calcHouseholderVectorCol(H, col + 1, col, H.numRows() - col - 1);
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

describe('Lower hessenberg zeroing', () => {
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
            let Q_k = makeHouseholderMatrix(v, H.numCols() - v.size(), H.numCols());
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

describe('Hessenbergization', () => {
});

describe('Tridiagonalization', () => {
    let A = new Matrix([
        4, 1, -2, 2,
        1, 2, 0, 1,
        -2, 0, 3, -2,
        2, 1, -2, -1], 4, 4);
    // test symmetric householder by generating hessenberg matrix QAQT = H
    test('Householder', () => {
        let Q = Matrix.identity(A.numRows());
        let T = A.clone();
        for (let iter = 0; iter + 2 < T.numCols(); ++iter) {
            let v = calcHouseholderVectorRow(T, iter, iter + 1);
            applyHouseholderFromLeft(v, Q, iter + 1);
            applyHouseholderFromLeft(v, T, iter + 1);
            applyHouseholderFromRight(v, T, iter + 1);
            // check matrices
        }
        expect(T.isHessenberg()).toBeTruthy();
        expect(T.isTridiagonal()).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(Matrix.mul(Q, Matrix.mul(A, Q.transpose())), T)).toBeLessThan(SmallTolerance);
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

    test('Hessenberg partial default', () => {
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

describe.skip('Triangular performance', () => {
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

describe('Hessenberg form decomposition', () => {
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
    const testData: Matrix[] = [
        new Matrix([
            1, 2, 3, 4,
            5, 6, 7, 8,
            3, 4, 2, -2,
            3, 5, 1, 2], 4, 4),
        new Matrix([
            0, -0.5, 0, 0,
            0, 0.25, -0.5, 0,
            0, -0.125, 0.25, -0.5,
            0, 0.0625, -0.125, 0.250], 4, 4),
        new Matrix([
            -0.5, 0, 0, 0,
            0.25, 0, -0.5, 0,
            -0.125, 0, 0.25, -0.5,
            0.0625, 0, -0.125, 0.250], 4, 4),
        new Matrix([
            0, -0.5, 0, 0,
            -0.5, 0.25, 0, 0,
            0.25, -0.125, 0, -0.5,
            -0.125, 0.0625, 0, 0.250], 4, 4),
        new Matrix([
            0, -0.5, 0, 0,
            0, 0.25, -0.5, 0,
            -0.5, -0.125, 0.25, 0,
            0.250, 0.0625, -0.125, 0], 4, 4)
    ];
    test.each(testData)('Hessenberg %#', (A: Matrix) => {
        let Q: Matrix = Matrix.empty(4, 4);
        let H = makeHessenberg(A, Q);
        expect(H.isHessenberg(true)).toBeTruthy();
        expect(Q.isOrthogonal()).toBeTruthy();
        expect(Matrix.lInfDistance(Matrix.mul(Matrix.mul(Q, A), Q.transpose()), H)).toBeLessThan(SmallTolerance);
    });
});

describe.only('Triangulization', () => {
    const testResult = (decomposition: OrthogonalDecomposition, isUpper: boolean, postMul: boolean, expectedMatrix: Matrix) => {
        expect(decomposition.Q).not.toBeNull();
        expect(decomposition.T).not.toBeNull();
        expect(decomposition.Q.isOrthogonal()).toBeTruthy();
        expect(decomposition.T.isTriangular(isUpper)).toBeTruthy();
        if (decomposition.Rank < Math.max(expectedMatrix.numRows(), expectedMatrix.numCols()))
            expect(decomposition.Q.determinant() * decomposition.T.determinant()).toBeCloseTo(decomposition.determinant());
        const actualMatrix = postMul ? Matrix.mul(decomposition.Q, decomposition.T) : Matrix.mul(decomposition.T, decomposition.Q);
        if (decomposition.P != null)
            decomposition.P.inverse().permuteInplace(actualMatrix);
        expect(Matrix.lInfDistance(actualMatrix, expectedMatrix)).toBeLessThan(Tolerance);
        expect(decomposition.Rank).toBeLessThanOrEqual(Math.min(expectedMatrix.numRows(), expectedMatrix.numCols()));
        if (isUpper) {
            if (postMul) {
                // QR
                if (decomposition.Rank < Math.min(expectedMatrix.numRows(), expectedMatrix.numCols())) {
                    for (let row = decomposition.Rank; row < decomposition.T.numRows(); ++row)
                        expect(decomposition.T.getRow(row).l2Norm()).toBeLessThan(SmallestTolerance);
                }
            } else {
                // RQ
                throw new Error("Not implemented");
            }

        } else {
            if (postMul) {
                // QL
                throw new Error("Not implemented");
            } else {
                // LQ
                throw new Error("Not implemented");
            }

        }
    };

    interface Method {
        method: ZeroingMethod, name: String
    };
    /*
    let fullRank = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    ], 4, 4);*/
    let fullRank = new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        8, 10, 11, 12,
        13, 13, 15, 16
    ], 4, 4);
    let nearSingular = new Matrix([
        1, 2, 1, 2,
        1, 2, 1, 3,
        1, 2, 2, 1,
        1, 2 + SmallestTolerance, 0, 0
    ], 4, 4);
    interface TestData {
        matrix: Matrix;
        rhs: Vector;
        solution: Vector;
        pseudoInverse: Matrix;
        rank: number;
        isMinNorm: number;
    }
    const testData = [
        // Square non-singular
        {
            matrix: new Matrix([
                1, 2, 3, 4,
                5, 6, 7, 8,
                -1, -2, 0, 2,
                2, 1, 4, 8
            ], 4, 4),
            solution: new Vector([-5, 2, 1, 1]), // exact solution
            rhs: new Vector([6, 2, 3, 4]),
            pseudoInverse: new Matrix([
                -7 / 6, 0.5, 1 / 3, 0,
                1, -1, -2, 1,
                -0.5, 1.5, 3, -2,
                5 / 12, -0.75, 4 / 3, 1
            ], 4, 4),
            rank: 4,
            isMinNorm: true
        },
        // Square singular inconsistent
        {
            matrix: new Matrix([
                1, 2, 3, 4,
                5, 6, 7, 8,
                -1, -2, 0, 2,
                2, 0, 4, 8
            ], 4, 4),
            solution: new Vector([0, 0, 0]), // least squares solution, todo: find correct solution
            rhs: new Vector([20, 70, 3, 46]),
            pseudoInverse: new Matrix([
                -5 / 6, 1 / 6, -1 / 3, 1 / 3,
                179 / 504, 17 / 504, 5 / 126, -61 / 252,
                55 / 252, 1 / 252, 4 / 63, -11 / 126,
                41 / 504, -13 / 504, 11 / 126, 17 / 252
            ], 4, 4),
            rank: 3,
            isMinNorm: false
        },
        // Square singular consistent
        {
            matrix: new Matrix([
                1, 2, 3, 4,
                5, 6, 7, 8,
                -1, -2, 0, 2,
                2, 0, 4, 8
            ], 4, 4),
            solution: new Vector([1, 2, 3, 4]),// min norm solution, todo: find correct solution
            rhs: new Vector([30, 70, 3, 46]),
            pseudoInverse: new Matrix([
                -5 / 6, 1 / 6, -1 / 3, 1 / 3,
                179 / 504, 17 / 504, 5 / 126, -61 / 252,
                55 / 252, 1 / 252, 4 / 63, -11 / 126,
                41 / 504, -13 / 504, 11 / 126, 17 / 252
            ], 4, 4),
            rank: 3,
            isMinNorm: true
        },
        // wide consistent
        {
            matrix: new Matrix([
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                13, 12, 13, 14, 15
            ], 3, 5),
            solution: new Vector([1, 2, 3, 4, 5]),
            rhs: new Vector([55, 130, 207]), // min norm solution
            pseudoInverse: new Matrix([
                0.5, - 1, 0.5,
                -1.06, 1.26, -0.5,
                -0.47, 0.62, -0.25,
                0.12, -0.2, 0,
                0.71, -0.66, 0.25
            ], 5, 3),
            rank: 3,
            isMinNorm: true
        },
        // wide inconsistent
        {
            matrix: new Matrix([
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                7, 9, 11, 13, 15
            ], 3, 5),
            solution: new Vector([0, 0, 0]), // least squares solution, todo: find correct solution
            rhs: new Vector([1, 2, 4]),
            pseudoInverse: new Matrix([
                22 / 75, 17 / 75, -1 / 15,
                -1 / 6, 2 / 15, -1 / 30,
                -1 / 25, 1 / 25, 0,
                13 / 150, -4 / 75, 1 / 30,
                16 / 75, -11 / 75, 1 / 15
            ], 5, 3),
            rank: 2,
            isMinNorm: false
        },
        // long consistent
        {
            matrix: new Matrix([
                1, 6, 7,
                2, 7, 9,
                3, 8, 11,
                0, 5, 5,
                1, 1, 2
            ], 5, 3),
            solution: new Vector([3, 2, 1]), // min norm solution
            rhs: new Vector([22, 29, 36, 15, 7]),
            pseudoInverse: new Matrix([
                -47 / 360, 11 / 360, 23 / 120, -7 / 24, 29 / 180,
                37 / 360, -1 / 360, -13 / 120, 5 / 24, -19 / 180,
                -1 / 36, 1 / 36, 1 / 12, -1 / 12, 1 / 18
            ], 3, 5),
            rank: 2,
            isMinNorm: true
        },
        // long inconsistent
        {
            matrix: new Matrix([
                1, 6, 7,
                2, 7, 9,
                3, 8, 11,
                4, 9, 13,
                5, 10, 15], 5, 3),
            solution: new Vector([-29 / 15, 64 / 15, 7 / 3]), // least squares solution, todo: check
            rhs: new Vector([44, 43, 52, 61, 70]),
            pseudoInverse: new Matrix([
                -22 / 75, -1 / 6, -1 / 25, 13 / 150, 16 / 75,
                17 / 75, 2 / 15, 1 / 25, -4 / 75, -11 / 75,
                -1 / 15, -1 / 30, 0, 1 / 30, 1 / 15
            ], 3, 5),
            rank: 2,
            isMinNorm: false
        },
        // wide consistent padded with zeroes
        {
            matrix: new Matrix([
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                13, 12, 13, 14, 15,
                0, 0, 0, 0, 0,
                0, 0, 0, 0, 0
            ], 5, 5),
            solution: new Vector([1, 2, 3, 4, 5]),
            rhs: new Vector([55, 130, 207, 0, 0]),
            pseudoInverse: new Matrix([], 5, 5),
            rank: 2,
            isMinNorm: true
        },
        // wide inconsistent padded with zeroes
        {
            matrix: new Matrix([
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                7, 9, 11, 13, 15,
                0, 0, 0, 0, 0,
                0, 0, 0, 0, 0
            ], 5, 5),
            solution: new Vector([]),
            rhs: new Vector([]),
            pseudoInverse: new Matrix([], 5, 5),
            rank: 2,
            isMinNorm: false
        },
        // long consistent padded with zeroes
        {
            matrix: new Matrix([
                1, 6, 7, 0, 0,
                2, 7, 9, 0, 0,
                3, 8, 11, 0, 0,
                0, 5, 5, 0, 0,
                1, 1, 2, 0, 0
            ], 5, 5),
            solution: new Vector([]),
            rhs: new Vector([]),
            pseudoInverse: new Matrix([
                -47 / 360, 11 / 360, 23 / 120, -7 / 24, 29 / 180,
                37 / 360, -1 / 360, -13 / 120, 5 / 24, -19 / 180,
                -1 / 36, 1 / 36, 1 / 12, -1 / 12, 1 / 18,
                0, 0, 0, 0, 0,
                0, 0, 0, 0, 0
            ], 5, 5),
            rank: 2,
            isMinNorm: true
        },
        // long inconsistent padded with zeroes
        {
            matrix: new Matrix([
                1, 6, 7, 0, 0,
                2, 7, 9, 0, 0,
                3, 8, 11, 0, 0,
                4, 9, 13, 0, 0,
                5, 10, 15, 0, 0
            ], 5, 5),
            solution: new Vector([]),
            rhs: new Vector([]),
            pseudoInverse: new Matrix([
                -22 / 75, -1 / 6, -1 / 25, 13 / 150, 16 / 75,
                17 / 75, 2 / 15, 1 / 25, -4 / 75, -11 / 75,
                -1 / 15, -1 / 30, 0, 1 / 30, 1 / 15,
                0, 0, 0, 0, 0,
                0, 0, 0, 0, 0
            ], 5, 5),
            rank: 2,
            isMinNorm: false
        }
    ];

    // todo:
    // Rank 2
    /*
    let rectWide = new Matrix([
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15
    ], 3, 5);*/
    // Rank 3
    let rectWide = new Matrix([
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        13, 12, 13, 14, 15
    ], 3, 5);
    let rectWideNearSingular = new Matrix([
        1, 2, 1, 2, 3,
        1, 2, 1, 3, 4,
        1, 2 + SmallestTolerance, 2, 1, 1
    ], 3, 5);
    let rectLong = rectWide.transpose();
    let rectLongNearSingular = rectWideNearSingular.transpose();
    test("Stuff", () => {
        let rectWideSqr = new Matrix([
            1, 2, 3, 4, 5,
            6, 7, 8, 9, 10,
            13, 12, 13, 14, 15,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
        ], 5, 5);
        let rhsWide = new Vector([1, 2, 3]);
        let rhsLong = new Vector([5, 4, 3, 4, 5]);
        let minSqrNormSolution = new Vector([0, -0.4 + 0.07 + 0.32, 0.6 - 0.14 - 0.48, 0.07, 0.16]);
        let leastSqrsSolution = new Vector([0.88, -1.47, 1]);
        let rectLongSqr = rectWideSqr.transpose();
        let decomposition = new QRTest(null, new OrthogonalDecompositionParams().setMethod(ZeroingMethod.Housholder).setCompact(false).setPivoting(true));
        decomposition.factorize(rectWideSqr);
        console.log(`rectWideSqr Q:${decomposition.Q.toString()}, R:${decomposition.R.toString()}, p:${decomposition.P.toString()}`);
        console.log(decomposition.P.inverse().permuteInplace(Matrix.mul(decomposition.Q, decomposition.R)).toString());

        decomposition.factorize(rectWide);
        console.log(`rectWide Q:${decomposition.Q.toString()}, R:${decomposition.R.toString()}, p:${decomposition.P.toString()}`);
        console.log(decomposition.P.inverse().permuteInplace(Matrix.mul(decomposition.Q, decomposition.R)).toString());
        /*
            decomposition.factorize(rectLongSqr);
            console.log(`rectLongSqr Q:${decomposition.Q.toString()}, R:${decomposition.R.toString()}, p:${decomposition.P.toString()}`);
            console.log(decomposition.P.inverse().permuteInplace(Matrix.mul(decomposition.Q, decomposition.R)).toString());
    
            decomposition.factorize(rectLong);
            console.log(`rectLong Q:${decomposition.Q.toString()}, R:${decomposition.R.toString()}, p:${decomposition.P.toString()}`);
            console.log(decomposition.P.inverse().permuteInplace(Matrix.mul(decomposition.Q, decomposition.R)).toString());
        */
    });

    describe.skip.each([{ method: ZeroingMethod.Givens, name: "Givens" }, { method: ZeroingMethod.Housholder, name: "Housholder" }])("Zeroing method %#",
        (method: Method) => {
            test('QR', () => {
                // square
                let withPivoting = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.QR, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(true).setCompact(false));
                testResult(withPivoting, true, true, fullRank);
                expect(withPivoting.Rank).toBe(4);

                let withoutPivotinig = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.QR, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(false).setCompact(false));
                testResult(withoutPivotinig, true, true, fullRank);
                expect(withoutPivotinig.Rank).toBe(4);

                withPivoting.factorize(nearSingular);
                testResult(withPivoting, true, true, nearSingular);
                expect(withPivoting.Rank).toBe(3);

                withoutPivotinig.factorize(nearSingular);
                expect(withoutPivotinig.Rank).toBe(0);
                expect(withoutPivotinig.Q).toBeNull();
                expect(withoutPivotinig.T).toBeNull();

                // rectangular
                /*withPivoting.factorize(rectLong);
    
    
                withPivoting.factorize(rectWide);
                withPivoting.factorize(rectLongNearSingular);
                withPivoting.factorize(rectWideNearSingular);*/
            });
            test.skip('RQ', () => {
                let withPivoting = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.RQ, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(true).setCompact(false));
                let withoutPivotinig = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.RQ, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(false).setCompact(false));
                testResult(withPivoting, true, false, fullRank);
                testResult(withoutPivotinig, true, false, fullRank);
                withPivoting.factorize(nearSingular);
                withoutPivotinig.factorize(nearSingular);
                testResult(withPivoting, true, false, nearSingular);
                expect(withoutPivotinig.Q).toBeNull();
                expect(withoutPivotinig.T).toBeNull();
            });
            test.skip("QL", () => {
                let withPivoting = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.QL, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(true).setCompact(false));
                let withoutPivotinig = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.QL, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(false).setCompact(false));
                testResult(withPivoting, false, true, fullRank);
                testResult(withoutPivotinig, false, true, fullRank);
                withPivoting.factorize(nearSingular);
                withoutPivotinig.factorize(nearSingular);
                testResult(withPivoting, false, true, nearSingular);
                expect(withoutPivotinig.Q).toBeNull();
                expect(withoutPivotinig.T).toBeNull();
            });
            test.skip("LQ", () => {
                let withPivoting = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.LQ, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(true).setCompact(false));
                let withoutPivotinig = new OrthogonalDecomposition(fullRank, OrthogonalDecompositionType.LQ, new OrthogonalDecompositionParams().setMethod(method.method).setPivoting(false).setCompact(false));
                testResult(withPivoting, false, false, fullRank);
                testResult(withoutPivotinig, false, false, fullRank);
                withPivoting.factorize(nearSingular);
                withoutPivotinig.factorize(nearSingular);
                testResult(withPivoting, false, false, nearSingular);
                expect(withoutPivotinig.Q).toBeNull();
                expect(withoutPivotinig.T).toBeNull();
            });
        });
    const data: Matrix[] = [new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        3, 4, 2, -2,
        3, 5, 1, 2], 4, 4),
    /*new Matrix([
        0, 2, 0, 4,
        0, 6, 0, 8,
        0, 4, 0, -2,
        0, 5, 0, 2], 4, 4),
    new Matrix([
        1, 2, 3, 4,
        5, 6, 7, 8,
        2, 4, 6, 8,
        5, 6, 7, 8], 4, 4),
    new Matrix([
        1, 2, 3, 4,
        0, 0, 0, 0,
        0, 0, 0, 0,
    5, 6, 7, 8], 4, 4)*/];
    describe.skip.each(data)("Square %#", (matrix: Matrix) => {
        describe.each([{ method: ZeroingMethod.Givens, name: "Givens" }, { method: ZeroingMethod.Housholder, name: "Housholder" }])("Zeroing method %#",
            (method: Method) => {
                test("QR", () => {
                    let decomposition = new OrthogonalDecomposition(matrix, OrthogonalDecompositionType.QR, new OrthogonalDecompositionParams().setMethod(method.method).setCompact(false));
                    testResult(decomposition, true, true, matrix);
                });
                test("RQ", () => {
                    let decomposition = new OrthogonalDecomposition(matrix, OrthogonalDecompositionType.RQ, new OrthogonalDecompositionParams().setMethod(method.method).setCompact(false));
                    testResult(decomposition, true, false, matrix);
                });
                test("QL", () => {
                    let decomposition = new OrthogonalDecomposition(matrix, OrthogonalDecompositionType.QL, new OrthogonalDecompositionParams().setMethod(method.method).setCompact(false));
                    testResult(decomposition, false, true, matrix);
                });
                test("LQ", () => {
                    let decomposition = new OrthogonalDecomposition(matrix, OrthogonalDecompositionType.LQ, new OrthogonalDecompositionParams().setMethod(method.method).setCompact(false));
                    testResult(decomposition, false, false, matrix);
                });
            });
    });
})