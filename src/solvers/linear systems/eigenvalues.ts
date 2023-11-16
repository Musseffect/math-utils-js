


// computes real port of eigenvalues using schur decomposition

import Matrix from "../../denseMatrix";
import { assert, sign } from "../../utils";
import Vector from "../../vector";
import { ConvergenseFailureException } from "./exceptions";

export interface givensCoeffs {
    c: number, s: number, r: number
};

export function givens(x: number, y: number): givensCoeffs {
    let c = 0.0;
    let s = 0.0;
    let r = 0.0;
    if (y == 0) {
        c = Math.sign(x);
        r = Math.abs(x);
    } else if (x == 0) {
        s = -Math.sign(y);
        r = Math.abs(y);
    } else if (Math.abs(x) > Math.abs(y)) {
        let t = y / x;
        let u = Math.sign(x) * Math.sqrt(1 + t * t);
        c = 1.0 / u;
        s = -c * t;
        r = x * u;
    } else {
        let t = x / y;
        let u = Math.sign(y) * Math.sqrt(1 + t * t);
        s = -1.0 / u;
        c = t / u;
        r = y * u;
    }
    return { c, s, r };
}
/** Multiply A by rotation matrix
 *        j ... i
 *  j |   c    -s
 *  . |
 *  i |   s     c
*/
export function applyGivensFromLeft(A: Matrix, givens: givensCoeffs, i: number, j: number): void {
    assert(i > j, "Incorrect order of indices");
    assert(A.numRows() > i, "Incorrect number of rows");
    for (let col = 0; col < A.numCols(); ++col) {
        let a = A.get(j, col);
        let b = A.get(i, col);
        A.set(j, col, givens.c * a - givens.s * b);
        A.set(i, col, givens.s * a + givens.c * b);
    }
}
/** Multiply A by transposed rotation matrix
 *        j ... i
*  j |    c     s
*  . |
*  i |   -s     c
*/
export function applyTransposeGivensFromLeft(A: Matrix, givens: givensCoeffs, i: number, j: number): void {
    assert(i > j, "Incorrect order of indices");
    assert(A.numRows() > i, "Incorrect number of rows");
    for (let col = 0; col < A.numCols(); ++col) {
        let a = A.get(j, col);
        let b = A.get(i, col);
        A.set(j, col, givens.c * a + givens.s * b);
        A.set(i, col, -givens.s * a + givens.c * b);
    }
}
/** Multiply A by transposed rotation matrix
 *        j ... i
*  j |    c     s
*  . |
*  i |   -s     c
*/
export function applyTransposeGivensFromRight(A: Matrix, givens: givensCoeffs, i: number, j: number) {
    assert(i > j, "Incorrect order of indices");
    assert(A.numCols() > i, "Incorrect number of cols");
    for (let row = 0; row < A.numRows(); ++row) {
        let a = A.get(row, j);
        let b = A.get(row, i);
        A.set(row, j, givens.c * a - givens.s * b);
        A.set(row, i, givens.s * a + givens.c * b);
    }
}
/** Multiply A by rotation matrix
 *        j ... i
 *  j |   c    -s
 *  . |
 *  i |   s     c
*/
export function applyGivensFromRight(A: Matrix, givens: givensCoeffs, i: number, j: number) {
    assert(i > j, "Incorrect order of indices");
    assert(A.numCols() > i, "Incorrect number of cols");
    for (let row = 0; row < A.numRows(); ++row) {
        let a = A.get(row, j);
        let b = A.get(row, i);
        A.set(row, j, givens.c * a + givens.s * b);
        A.set(row, i, -givens.s * a + givens.c * b);
    }
}
/**
 *  Make matrix
 * 
 * [c, -s]\
 * [s, c]
 *  */

export function makeGivensMatrix(givens: givensCoeffs, matrixSize: number, i: number, j: number): Matrix {
    assert(i > j, "Incorrect order of indices");
    assert(matrixSize > i, "Incorrect matrix size");
    let m = Matrix.empty(matrixSize, matrixSize);
    m.set(j, j, givens.c);
    m.set(j, i, -givens.s);
    m.set(i, i, givens.c);
    m.set(i, j, givens.s);
    for (let row = 0; row < matrixSize; ++row) {
        if (row != i && row != j)
            m.set(row, row, 1);
    }
    return m;
}

function processHouseholderVectorInplace(v: Vector): Vector {
    let ro = -sign(v.get(0));
    let xNormSqr = v.squaredLength();
    let xNorm = Math.sqrt(xNormSqr);
    let firstElement = v.get(0);
    v.set(0, v.get(0) - ro * xNorm);
    v.scaleSelf(1.0 / Math.sqrt(xNormSqr - firstElement * firstElement + v.get(0) * v.get(0)));
    return v;
}

export function calcHouseholderVectorCol(A: Matrix, row: number, col: number): Vector {
    assert(row < A.numRows(), "Incorrect row");
    let v = A.subColumn(row, col, A.numRows() - row);
    return processHouseholderVectorInplace(v);
}

export function calcHouseholderVectorRow(A: Matrix, row: number, col: number): Vector {
    assert(col < A.numCols(), "Incorrect row");
    let v = A.subRow(row, col, A.numCols() - col);
    return processHouseholderVectorInplace(v);
}

// should be 4/3N^3 + O(N^2)
export function makeTridiagonalInplace(A: Matrix, Q?: Matrix): Matrix {
    assert(A.isSymmetric(), "A must be symmetric");
    if (Q) {
        Q.setFromMatrix(Matrix.identity(A.numRows()));
    }
    for (let outerCol = 0; outerCol + 2 < A.numCols(); ++outerCol) {
        let shift = outerCol + 1;
        let u = A.subRow(outerCol, shift, A.numRows() - shift);
        let xNormSqr = u.squaredLength();
        let xNorm = Math.sqrt(xNormSqr);
        let firstElement = u.get(0);
        let ro = -sign(firstElement);
        // first element of the column is alpha, elements below are zero 
        let alpha = ro * xNorm;
        u.set(0, firstElement - alpha);
        u.scaleSelf(1.0 / Math.sqrt(xNormSqr - firstElement * firstElement + u.get(0) * u.get(0)));
        // premultiply all rows
        // first (col + 1) rows won't change
        // set first column
        A.set(shift, outerCol, alpha);
        A.set(outerCol, shift, alpha);
        for (let row = shift + 1; row < A.numRows(); ++row) {
            A.set(row, outerCol, 0);
            A.set(outerCol, row, 0);
        }
        // calc pre and post multiplicaiton simultaniously
        let v = Vector.empty(A.numRows() - shift);
        // calc 2 * A * u and store in v
        for (let row = shift; row < A.numRows(); ++row) {
            let value = 0.0;
            for (let col = shift; col < A.numCols(); ++col) {
                value += A.get(row, col) * u.get(col - shift);
            }
            v.set(row - shift, 2 * value);
        }
        let utAu = Vector.dot(u, v);
        // calc v = 2 * A * u - 2u * (ut * A * u)
        for (let idx = 0; idx < v.size(); ++idx) {
            v.set(idx, v.get(idx) - u.get(idx) * utAu);
        }

        // calc A - uvT - vutv
        for (let row = shift; row < A.numRows(); ++row) {
            for (let col = row; col < A.numCols(); ++col) {
                A.set(row, col, A.get(row, col) - v.get(row - shift) * u.get(col - shift) - v.get(col - shift) * u.get(row - shift));
            }
        }
        // preserve symmetry
        for (let row = shift; row < A.numRows(); ++row) {
            for (let col = row + 1; col < A.numCols(); ++col)
                A.set(col, row, A.get(row, col));
        }
        if (Q) {
            applyHouseholderFromLeft(u, Q, shift);
        }
    }
    return A;
}

export function makeTridiagonalInplaceAlt(A: Matrix, Q?: Matrix): Matrix {
    assert(A.isSymmetric(), "A must be symmetric");
    for (let outerCol = 0; outerCol + 2 < A.numCols(); ++outerCol) {
        let shift = outerCol + 1;
        let xNormSqr = 0;
        for (let col = shift; col < A.numCols(); ++col)
            xNormSqr += Math.pow(A.get(outerCol, col), 2);
        let xNorm = Math.sqrt(xNormSqr);
        let firstElement = A.get(outerCol, shift);
        let ro = -sign(firstElement);
        // first element of the column is alpha, elements below are zero 
        let alpha = ro * xNorm;
        let newFirstElement = firstElement - alpha;
        A.set(outerCol, shift, firstElement - alpha);
        let newLength = Math.sqrt(xNormSqr - firstElement * firstElement + newFirstElement * newFirstElement);
        for (let col = shift; col < A.numCols(); ++col)
            A.set(outerCol, col, A.get(outerCol, col) / newLength);
        //let u = A.subRow(outerCol, shift, A.numRows() - shift);
        // premultiply all rows
        // first (col + 1) rows won't change
        // set first column
        A.set(shift, outerCol, alpha);
        for (let row = shift + 1; row < A.numRows(); ++row)
            A.set(row, outerCol, 0);
        // calc pre and post multiplicaiton simultaniously
        let v = Vector.empty(A.numRows() - shift);
        // calc 2 * A * u and store in v
        for (let row = shift; row < A.numRows(); ++row) {
            let value = 0.0;
            for (let col = shift; col < A.numCols(); ++col) {
                value += A.get(row, col) * /*u*/A.get(outerCol, col);
            }
            v.set(row - shift, 2 * value);
        }
        let utAu = 0;
        for (let col = shift; col < A.numCols(); ++col)
            utAu += /*u*/A.get(outerCol, col) * v.get(col - shift);
        // calc v = 2 * A * u - 2u * (ut * A * u)
        for (let idx = 0; idx < v.size(); ++idx) {
            v.set(idx, v.get(idx) - /*u*/A.get(outerCol, idx + shift) * utAu);
        }

        // calc A - uvT - vutv
        for (let row = shift; row < A.numRows(); ++row) {
            for (let col = row; col < A.numCols(); ++col) {
                A.set(row, col, A.get(row, col) - v.get(row - shift) * /*u*/A.get(outerCol, col) - v.get(col - shift) * /*u*/A.get(outerCol, row));
            }
        }
        // preserve symmetry
        for (let row = shift; row < A.numRows(); ++row) {
            for (let col = row + 1; col < A.numCols(); ++col)
                A.set(col, row, A.get(row, col));
        }
    }
    if (Q) {
        // Set initial Q with vlaues for 2x2 Q_k
        let v1 = A.get(A.numRows() - 3, A.numCols() - 2);
        let v2 = A.get(A.numRows() - 3, A.numCols() - 1);
        let v22 = v2 * v2;
        let v11 = v1 * v1;
        let v12 = v1 * v2;
        Q.set(A.numRows() - 2, A.numCols() - 2, 1 - 2 * v11);
        Q.set(A.numRows() - 1, A.numCols() - 2, - 2 * v12);
        Q.set(A.numRows() - 2, A.numCols() - 1, - 2 * v12);
        Q.set(A.numRows() - 1, A.numCols() - 1, 1 - 2 * v22);

        Q.set(A.numRows() - 3, A.numCols() - 3, 1);
        // accumulate v in A and construct Q at the end by multiplying from last to first
        for (let aRow = A.numCols() - 4; aRow >= 0; --aRow) {
            let qCol = aRow + 1;
            Q.set(aRow, aRow, 1);
            let v1_2 = 2 * A.get(aRow, qCol);
            // set the first row of non-identity part of Q
            for (let col = qCol; col < A.numCols(); ++col) {
                Q.set(qCol, col, Q.get(qCol, col) - v1_2 * A.get(aRow, col));
            }
            for (let row = qCol + 1; row < A.numRows(); ++row) {
                let vDotX = 0.0;
                for (let col = qCol + 1; col < A.numCols(); ++col) {
                    vDotX += A.get(aRow, col) * Q.get(row, col);
                }
                vDotX *= 2;
                for (let col = qCol; col < A.numCols(); ++col) {
                    Q.set(row, col, Q.get(row, col) - A.get(aRow, col) * vDotX);
                }
            }
        }
    }
    // set correct values for upper triangle of A
    for (let row = 0; row + 2 < A.numRows(); ++row)
        for (let col = row + 1; col < A.numCols(); ++col)
            A.set(row, col, A.get(col, row));
    return A;
}

export function makeTridiagonal(A: Matrix, Q?: Matrix) {
    return makeTridiagonalInplace(A.clone(), Q,);
}

export function makeTridiagonalAlt(A: Matrix, Q?: Matrix) {
    return makeTridiagonalInplaceAlt(A.clone(), Q,);
}

// todo: test
export function applyHouseholderFromLeft(v: Vector, A: Matrix, idx: number) {
    for (let col = 0; col < A.numCols(); ++col) {
        let vDotX = 0.0;
        for (let row = idx; row < A.numRows(); ++row) {
            vDotX += v.get(row - idx) * A.get(row, col);
        }
        vDotX *= 2;
        for (let row = idx; row < A.numRows(); ++row) {
            A.set(row, col, A.get(row, col) - v.get(row - idx) * vDotX);
        }
    }
}

// todo: test
export function applyHouseholderFromRight(v: Vector, A: Matrix, idx: number) {
    for (let row = 0; row < A.numRows(); ++row) {
        let vDotX = 0.0;
        for (let col = idx; col < A.numCols(); ++col)
            vDotX += v.get(col - idx) * A.get(row, col);
        vDotX *= 2;
        for (let col = idx; col < A.numCols(); ++col)
            A.set(row, col, A.get(row, col) - v.get(col - idx) * vDotX);
    }
}

export function makeHouseholder(v: Vector, size: number): Matrix {
    let result = Matrix.identity(size);
    let idx = size - v.size();
    for (let row = 0; row < v.size(); ++row) {
        for (let col = 0; col < v.size(); ++col)
            result.set(row + idx, col + idx, result.get(row + idx, col + idx) - 2 * v.get(row) * v.get(col));
    }
    return result;
}

// todo: test, should be 6n^3 + O(n^2)
export function makeHessenbergInplace(A: Matrix, Q?: Matrix): Matrix {
    //if (A.isSymmetric()) return makeTridiagonalInplace(A, Q);
    assert(A.isSquare(), "Non-square matrix");
    // todo: store n-1 elements of v in H and reconstruct Q afterwards: |v| = 1 by applying householder from the right
    // Q = P1*P2...PN-2
    if (Q) {
        Q.setFromMatrix(Matrix.identity(A.numRows()));
    }
    for (let outerCol = 0; outerCol + 2 < A.numCols(); ++outerCol) {
        let shift = outerCol + 1;
        let v = A.subColumn(shift, outerCol, A.numRows() - shift);
        let xNormSqr = v.squaredLength();
        let xNorm = Math.sqrt(xNormSqr);
        let firstElement = v.get(0);
        let ro = -sign(firstElement);
        // first element of the column is alpha, elements below are zero 
        let alpha = ro * xNorm;
        v.set(0, firstElement - alpha);
        v.scaleSelf(1.0 / Math.sqrt(xNormSqr - firstElement * firstElement + v.get(0) * v.get(0)));
        // premultiply all rows
        // first (col + 1) rows won't change
        // set first column
        A.set(shift, outerCol, alpha);
        for (let row = shift + 1; row < A.numRows(); ++row)
            A.set(row, outerCol, 0);
        // set other columns ~O(2N^2)
        for (let col = shift; col < A.numCols(); ++col) {
            let vDotX = 0.0;
            for (let row = shift; row < A.numRows(); ++row) {
                vDotX += v.get(row - shift) * A.get(row, col);
            }
            vDotX *= 2;
            for (let row = shift; row < A.numRows(); ++row) {
                A.set(row, col, A.get(row, col) - v.get(row - shift) * vDotX);
            }
        }
        for (let row = 0; row < A.numRows(); ++row) {
            let vDotX = 0.0;
            for (let col = shift; col < A.numCols(); ++col) {
                vDotX += v.get(col - shift) * A.get(row, col);
            }
            vDotX *= 2;
            for (let col = shift; col < A.numCols(); ++col) {
                A.set(row, col, A.get(row, col) - v.get(col - shift) * vDotX);
            }
        }
        // todo: accumulate v in A and construct Q at the end by multiplying from last to first;
        if (Q) {
            applyHouseholderFromLeft(v, Q, shift);
        }
    }
    return A;
}

/**
 * Compute QHQT = A
 */
export function makeHessenberg(A: Matrix, Q?: Matrix): Matrix {
    return makeHessenbergInplace(A.clone(), Q);
}

function makeHessenbergInplaceAltWithQ(A: Matrix, Q: Matrix): Matrix {
    Q.setFromMatrix(Matrix.identity(A.numRows()));
    if (A.numCols() < 3) return A;
    for (let outerCol = 0; outerCol + 2 < A.numCols(); ++outerCol) {
        let shift = outerCol + 1;
        for (let row = shift; row < A.numRows(); ++row)
            Q.set(shift, row, A.get(row, outerCol));
        let xNormSqr = 0.0;
        for (let row = shift; row < A.numRows(); ++row)
            xNormSqr += Math.pow(A.get(row, outerCol), 2);
        let xNorm = Math.sqrt(xNormSqr);
        let firstElement = Q.get(shift, shift);
        let ro = -sign(firstElement);
        // first element of the column is alpha, elements below are zero 
        let alpha = ro * xNorm;
        let newFirstElement = firstElement - alpha;
        Q.set(shift, shift, newFirstElement);
        const factor = 1.0 / Math.sqrt(xNormSqr - firstElement * firstElement + newFirstElement * newFirstElement);
        for (let row = shift; row < A.numRows(); ++row)
            Q.set(shift, row, /*u*/Q.get(shift, row) * factor);
        // premultiply all rows
        // first (col + 1) rows won't change
        // set first column
        A.set(shift, outerCol, alpha);
        for (let row = shift + 1; row < A.numRows(); ++row)
            A.set(row, outerCol, 0);
        // set other columns ~O(2N^2)
        for (let col = shift; col < A.numCols(); ++col) {
            let vDotX = 0.0;
            for (let row = shift; row < A.numRows(); ++row) {
                vDotX += /*u*/Q.get(shift, row) * A.get(row, col);
            }
            vDotX *= 2;
            for (let row = shift; row < A.numRows(); ++row) {
                A.set(row, col, A.get(row, col) - /*u*/Q.get(shift, row) * vDotX);
            }
        }
        for (let row = 0; row < A.numRows(); ++row) {
            let vDotX = 0.0;
            for (let col = shift; col < A.numCols(); ++col) {
                vDotX += /*u*/Q.get(shift, col) * A.get(row, col);
            }
            vDotX *= 2;
            for (let col = shift; col < A.numCols(); ++col) {
                A.set(row, col, A.get(row, col) - /*u*/Q.get(shift, col) * vDotX);
            }
        }
    }
    // set first 2x2 non identity part of matrix Q
    let v1 = Q.get(Q.numRows() - 2, Q.numCols() - 2);
    let v2 = Q.get(Q.numRows() - 2, Q.numCols() - 1);
    let v22 = v2 * v2;
    let v11 = v1 * v1;
    let v12 = v1 * v2;
    Q.set(Q.numRows() - 2, Q.numCols() - 2, 1 - 2 * v11);
    Q.set(Q.numRows() - 1, Q.numCols() - 2, - 2 * v12);
    Q.set(Q.numRows() - 2, Q.numCols() - 1, - 2 * v12);
    Q.set(Q.numRows() - 1, Q.numCols() - 1, 1 - 2 * v22);
    // accumulate v in A and construct Q at the end by multiplying from last to first
    for (let aCol = Q.numCols() - 4; aCol >= 0; --aCol) {
        let qCol = aCol + 1;
        let v1 = Q.get(qCol, qCol);
        let v1_2 = 2 * v1;
        for (let row = qCol + 1; row < Q.numRows(); ++row) {
            let vDotX = 0.0;
            for (let col = qCol + 1; col < Q.numCols(); ++col) {
                vDotX += /*u*/Q.get(qCol, col) * Q.get(row, col);
            }
            vDotX *= 2;
            for (let col = qCol; col < Q.numCols(); ++col) {
                Q.set(row, col, Q.get(row, col) - /*u*/Q.get(qCol, col) * vDotX);
            }
        }
        Q.set(qCol, qCol, 1 - v1 * v1_2);
        // set the first row of non-identity part of Q
        for (let col = qCol + 1; col < Q.numCols(); ++col) {
            Q.set(qCol, col, - v1_2 * /*u*/Q.get(qCol, col));
        }
    }
    return A;
}

// todo: test
export function makeHessenbergInplaceAlt(A: Matrix, Q?: Matrix): Matrix {
    //if (A.isSymmetric()) return makeTridiagonalInplace(A, Q);
    assert(A.isSquare(), "Non-square matrix");
    // todo: store n-1 elements of v in H and reconstruct Q afterwards: |v| = 1 by applying householder from the right
    // Q = P1*P2...PN-2
    if (Q) return makeHessenbergInplaceAltWithQ(A, Q);
    for (let outerCol = 0; outerCol + 2 < A.numCols(); ++outerCol) {
        let shift = outerCol + 1;
        let v = A.subColumn(shift, outerCol, A.numRows() - shift);
        let xNormSqr = v.squaredLength();
        let xNorm = Math.sqrt(xNormSqr);
        let firstElement = v.get(0);
        let ro = -sign(firstElement);
        // first element of the column is alpha, elements below are zero 
        let alpha = ro * xNorm;
        v.set(0, firstElement - alpha);
        v.scaleSelf(1.0 / Math.sqrt(xNormSqr - firstElement * firstElement + v.get(0) * v.get(0)));
        // premultiply all rows
        // first (col + 1) rows won't change
        // set first column
        A.set(shift, outerCol, alpha);
        for (let row = shift + 1; row < A.numRows(); ++row)
            A.set(row, outerCol, 0);
        // set other columns ~O(2N^2)
        for (let col = shift; col < A.numCols(); ++col) {
            let vDotX = 0.0;
            for (let row = shift; row < A.numRows(); ++row) {
                vDotX += v.get(row - shift) * A.get(row, col);
            }
            vDotX *= 2;
            for (let row = shift; row < A.numRows(); ++row) {
                A.set(row, col, A.get(row, col) - v.get(row - shift) * vDotX);
            }
        }
        for (let row = 0; row < A.numRows(); ++row) {
            let vDotX = 0.0;
            for (let col = shift; col < A.numCols(); ++col) {
                vDotX += v.get(col - shift) * A.get(row, col);
            }
            vDotX *= 2;
            for (let col = shift; col < A.numCols(); ++col) {
                A.set(row, col, A.get(row, col) - v.get(col - shift) * vDotX);
            }
        }
    }
    return A;
}

export function makeHessenbergAlt(A: Matrix, Q?: Matrix): Matrix {
    return makeHessenbergInplaceAlt(A.clone(), Q);
}
// todo: Givens rotations, complex eigenvalues
/**
 * Compute eigenvalues with Francis double step QR algorithm for arbitary square matrix
 * @param A Input matrix
 * @param numIters number of QR iterations
 * @param tolerance tolerance for zeroed elements
 * @returns array of real eigenvalues
 */
export function calcEigenvalues(A: Matrix, numIters: number, tolerance: number): number[] {
    assert(A.isSquare(), "Expected square matrix");
    let eigenvalues: number[] = new Array(A.numCols());
    if (A.numCols() == 1) {
        eigenvalues[0] = A.get(0, 0);
        return eigenvalues;
    }
    else if (A.numCols() == 2) {
        let b = A.get(0, 0) + A.get(1, 1);
        let c = A.get(0, 0) * A.get(1, 1) - A.get(0, 1) * A.get(1, 0);
        let D = (b * b - 4.0 * c);
        eigenvalues[0] = b * 0.5;
        eigenvalues[1] = b * 0.5;
        if (D > 0) {
            D = Math.sqrt(D);
            eigenvalues[0] += D * 0.5;
            eigenvalues[1] -= D * 0.5;
        }
        return eigenvalues;
    }
    A = A.clone();

    let u = Vector.empty(A.numCols());
    if (!A.isHessenberg()) {
        // turn into hessenberg matrix
        for (let i = 0; i < u.size() - 2; i++) {
            let xNorm = 0.0;
            for (let k = 0, j = i + 1; j < u.size(); j++, k++) {
                u.set(k, A.get(j, i));
                xNorm += u.get(k) * u.get(k);
            }
            let ro = -Math.sign(A.get(i + 1, i));
            let uNorm = xNorm - A.get(i + 1, i) * A.get(i + 1, i);
            u.set(0, u.get(0) - ro * Math.sqrt(xNorm));
            uNorm += u.get(0) * u.get(0);
            uNorm = Math.sqrt(uNorm);
            u.scaleSelf(1.0 / uNorm);
            let u_a = Vector.empty(u.size() - i); //uk* Ak+1:n,k:n
            // premultiply by Q_i
            for (let j = i; j < u.size(); j++) {
                let value = 0.0;
                for (let k = i + 1; k < u.size(); k++)
                    value += u.get(k - i - 1) * A.get(k, j);
                u_a.set(j - i, value);
            }

            for (let j = i + 1; j < u.size(); j++) {
                for (let k = i; k < u.size(); k++)
                    A.set(j, k, A.get(j, k) - u.get(j - i - 1) * 2.0 * u_a.get(k - i));
            }
            // postmultiply by Q_i
            u_a = Vector.empty(u.size());
            for (let j = 0; j < u.size(); j++) {
                let value = 0.0;
                for (let k = i + 1; k < u.size(); k++)
                    value += u.get(k - i - 1) * A.get(j, k);
                u_a.set(j, value);
            }

            for (let j = 0; j < u.size(); j++) {
                for (let k = i + 1; k < u.size(); k++)
                    A.set(j, k, A.get(j, k) - 2.0 * u_a.get(j) * u.get(k - i - 1));
            }
        }
    }
    //console.log(`Hessenberg ${A.toString()}`);
    // A = makeHessenberg(A);

    for (let i = 0; i < u.size() - 2; i++) {
        for (let j = i + 2; j < u.size(); j++)
            A.set(j, i, 0.0);
    }
    // Francis double step QR
    let iter = 0;
    for (let p = u.size() - 1; p > 1;) {
        let q = p - 1;
        let s = A.get(q, q) + A.get(p, p);
        let t = A.get(q, q) * A.get(p, p) - A.get(p, q) * A.get(q, p);
        let x = A.get(0, 0) * A.get(0, 0) + A.get(0, 1) * A.get(1, 0) - s * A.get(0, 0) + t;
        let y = A.get(1, 0) * (A.get(0, 0) + A.get(1, 1) - s);
        let z = A.get(1, 0) * A.get(2, 1);
        for (let k = 0; k <= p - 2; k++) {
            let r = Math.max(0, k - 1);
            let p_v = new Vector([x, y, z]);
            let ro = -Math.sign(x);
            p_v.set(0, p_v.get(0) - ro * p_v.l2Norm());
            p_v.normalize();

            let p_t = new Array(u.size() - r);
            for (let j = r, m = 0; j < u.size(); j++, m++) {
                let temp = 0.0;
                for (let i = k, l = 0; l < 3; i++, l++)
                    temp += p_v.get(l) * A.get(i, j);
                p_t[m] = temp;
            }
            for (let j = k, l = 0; l < 3; j++, l++) {
                for (let i = r; i < u.size(); i++)
                    A.set(j, i, A.get(j, i) - 2.0 * p_v.get(l) * p_t[i - r]);
            }
            r = Math.min(k + 3, p);
            p_t = new Array(r + 1);
            for (let j = 0; j <= r; j++) {
                let value = 0.0;
                for (let i = k, l = 0; l < 3; i++, l++)
                    value += p_v.get(l) * A.get(j, i);
                p_t[j] = value;
            }

            for (let i = 0; i <= r; i++) {
                for (let j = k, l = 0; l < 3; j++, l++)
                    A.set(i, j, A.get(i, j) - 2.0 * p_v.get(l) * p_t[i]);
            }
            x = A.get(k + 1, k);
            y = A.get(k + 2, k);
            if (k < p - 2)
                z = A.get(k + 3, k);
        }

        let p_v = new Vector([x, y]);
        let ro = -Math.sign(x);
        p_v.set(0, p_v.get(0) - ro * p_v.l2Norm());
        p_v.normalize();

        let p_t = new Array(u.size() - p + 2);
        for (let j = p - 2, m = 0; j < u.size(); j++, m++) {
            let temp = 0.0;
            for (let i = q; i <= p; i++)
                temp += p_v.get(i - q) * A.get(i, j);
            p_t[m] = temp;
        }
        for (let i = q; i <= p; i++) {
            for (let j = p - 2, m = 0; j < u.size(); j++, m++)
                A.set(i, j, A.get(i, j) - 2.0 * p_v.get(i - q) * p_t[m]);
        }


        p_t = new Array(p + 1);
        for (let j = 0; j <= p; j++) {
            let value = 0.0;
            for (let i = p - 1, m = 0; i <= p; i++, m++)
                value += p_v.get(m) * A.get(j, i);
            p_t[j] = value;
        }

        for (let i = 0; i <= p; i++) {
            for (let j = p - 1, m = 0; j <= p; j++, m++)
                A.set(i, j, A.get(i, j) - 2.0 * p_v.get(m) * p_t[i]);
        }
        if (Math.abs(A.get(p, q)) < tolerance * (Math.abs(A.get(q, q)) + Math.abs(A.get(p, p)))) {
            A.set(p, q, 0);
            p = p - 1;
            q = p - 1;
        } else if (Math.abs(A.get(p - 1, q - 1)) < tolerance * (Math.abs(A.get(q - 1, q - 1)) + Math.abs(A.get(q, q)))) {
            A.set(p - 1, q - 1, 0);
            p = p - 2;
            q = p - 1;
        }
        iter++;
        if (iter > numIters)
            throw new ConvergenseFailureException("EignevaluesSolver");
    }
    for (let i = 0; i < A.numCols(); i++) {
        if (i > 0 && Math.abs(A.get(i, i - 1)) > tolerance * 10.0) //complex eigenvalues
        {
            let b = A.get(i - 1, i - 1) + A.get(i, i);
            let c = A.get(i - 1, i - 1) * A.get(i, i) - A.get(i - 1, i) * A.get(i, i - 1);
            let D = b * b - 4.0 * c;
            let x1 = b * 0.5;
            let x2 = b * 0.5;
            if (D > 0) {
                D = Math.sqrt(D);
                x1 += D * 0.5;
                x2 -= D * 0.5;
            }
            eigenvalues[i - 1] = x1;
            eigenvalues[i] = x2;
        } else {
            eigenvalues[i] = A.get(i, i);
        }
    }
    return eigenvalues;
}

// todo: use for regularization of Newton optimization method by clamping negative eigenvalues
/**
 * Q - orthogonal matrix, D - diagonal matrix of eigenvalues of A
 */

class SymmetricEigendecomposition {
    private q: Matrix = null;
    private d: Vector = null;
    private A: Matrix = null;
    constructor(A: Matrix | null = null) {
        if (A == null) return;
        this.factorize(A);
    }
    public factorize(A: Matrix) {
        assert(A.isSymmetric(), "Expected symmetric matrix");
        this.A = A;
        this.q = null;
        this.d = null;
        let T = A.clone();
        this.q = Matrix.identity(A.numCols());
        if (!T.isTridiagonal())
            makeTridiagonalInplace(T, this.q);
        // symmetric matrices have only real eigenvalues so single shift algorithm can be used
        throw new Error("Not implemented");
    }
    public get D(): Vector {
        return this.d;
    }
    public get Q(): Matrix {
        return this.q;
    }
}

class Eigendecomposition {
    private q: Matrix = null;
    private d: Matrix = null;
    private A: Matrix = null;
    constructor(A: Matrix | null = null) {
        if (A == null) return;
        this.factorize(A);
    }
    public factorize(A: Matrix) {
        this.A = A;
        this.q = null;
        this.d = A.clone();
        this.q = Matrix.identity(A.numCols());
        if (!this.d.isHessenberg(true))
            makeHessenbergInplace(this.d, this.q);
        // symmetric matrices have only real eigenvalues so single shift algorithm can be used
        throw new Error("Not implemented");
    }
    public get D(): Matrix {
        return this.d;
    }
    public get Q(): Matrix {
        return this.q;
    }
}
export function calcEigendecomposition(A: Matrix): { Q: Matrix, D: Matrix } {
    throw new Error("Not implemented");
}

// todo: Jacobi method for eigenvalues of symmetric matrix, also for singular values too