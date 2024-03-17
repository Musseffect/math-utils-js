import Matrix from "../../denseMatrix";
import { assert, sign } from "../../utils";
import Vector from "../../vector";

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