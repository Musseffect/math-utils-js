//todo:

import Matrix from "../../denseMatrix";
import { DiagonalType, TriMatrixType, TriMatrixView } from "../../triMatrixView";
import { assert, near, SmallTolerance } from "../../utils";
import Vector from "../../vector";
import { NotPositiveDefiniteMatrixException } from "./exceptions";

const SolverName = "'Cholesky'";

export default class LLT {
    llt: Matrix;
    constructor() {
    }
    factorize(A: Matrix): void {
        assert(A.isSquare(), "Non-square matrix");
        this.llt = A.clone();
        let rank = A.width();
        for (let row = 0; row < rank; ++row) {
            for (let column = 0; column <= row; ++column) {
                let value = this.llt.get(row, column);
                for (let k = 0; k < column; ++k)
                    value -= this.llt.get(row, k) * this.llt.get(column, k);
                if (row == column) {
                    if (value < 0.0) throw new NotPositiveDefiniteMatrixException(SolverName);
                    value = Math.sqrt(value);
                }
                else {
                    value = value / this.llt.get(column, column);
                }
                this.llt.set(row, column, value);
                this.llt.set(column, row, value);
            }
        }
    };
    L(): TriMatrixView {
        return new TriMatrixView(this.llt, TriMatrixType.lower, DiagonalType.Unit);
    }
    LT(): TriMatrixView {
        return new TriMatrixView(this.llt, TriMatrixType.upper, DiagonalType.Existing);
    }
    solve(rhs: Vector): Vector {
        assert(rhs.size() == this.llt.width(), "Incompatible RHS");
        const rank = this.llt.width();
        let y = Vector.empty(this.llt.width());
        for (let row = 0; row < rank; ++row) {
            let value = rhs.get(row);
            for (let column = 0; column < row; ++column)
                value -= this.llt.get(row, column) * y.get(column);
            value /= this.llt.get(row, row);
            y.set(row, value);
        }
        let x = Vector.empty(rank);
        for (let row = rank - 1; row >= 0; --row) {
            let value = y.get(row);
            for (let column = row + 1; column < rank; ++column)
                value -= this.llt.get(row, column) * x.get(column);
            value /= this.llt.get(row, row);
            x.set(row, value);
        }
        return x;
    }
    solveMatrix(rhs: Matrix): Matrix {
        assert(rhs.height() == this.llt.width(), "Incompatible RHS");
        throw new Error("Not implemented");
    }
    inverse(): Matrix {
        throw new Error("Not implemented");
    }
    // Solve inplace
    static solve(A: Matrix, b: Vector) {
        assert(A.width() == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.isSquare(), "Non-square matrix");

        let rank = b.size();
        for (let row = 0; row < rank; ++row) {
            for (let column = 0; column <= row; ++column) {
                let value = A.get(row, column);
                for (let k = 0; k < column; ++k)
                    value -= A.get(row, k) * A.get(column, k);
                if (row == column) {
                    if (value < 0.0) throw new NotPositiveDefiniteMatrixException(SolverName);
                    value = Math.sqrt(value);
                }
                else {
                    value = value / A.get(column, column);
                }
                A.set(row, column, value);
                A.set(column, row, value);
            }
        }
        let y = Vector.empty(rank);
        for (let row = 0; row < rank; ++row) {
            let value = b.get(row);
            for (let column = 0; column < row; ++column)
                value -= A.get(row, column) * y.get(column);
            value /= A.get(row, row);
            y.set(row, value);
        }
        let x = Vector.empty(rank);
        for (let row = rank - 1; row >= 0; --row) {
            let value = y.get(row);
            for (let column = row + 1; column < rank; ++column)
                value -= A.get(row, column) * x.get(column);
            value /= A.get(row, row);
            x.set(row, value);
        }
        return x;
    }
}