//todo:

import Matrix from "../../denseMatrix";
import { TriMatrixType, TriMatrixView } from "../../triMatrixView";
import { assert, near, SmallTolerance } from "../../utils";
import Vector from "../../vector";
import { NotPositiveDefiniteMatrixException } from "./exceptions";

const SolverName = "'Cholesky'";

export default class LL {
    ll: Matrix;
    A: Matrix;
    private compute() {
        assert(this.A.isSquare(), "Non-square matrix");

        let rank = this.A.width();
        for (let row = 0; row < rank; ++row) {
            for (let column = 0; column <= row; ++column) {
                let value = this.A.get(row, column);
                for (let k = 0; k < column; ++k)
                    value -= this.ll.get(row, k) * this.ll.get(column, k);
                if (row == column) {
                    if (value < 0.0) throw new NotPositiveDefiniteMatrixException(SolverName);
                    value = Math.sqrt(value);
                }
                else {
                    value = value / this.ll.get(column, column);
                }
                this.ll.set(row, column, value);
                this.ll.set(column, row, value);
            }
        }
    }
    constructor(A: Matrix) {
        this.A = A;
        this.compute();
    }
    L(): TriMatrixView {
        return new TriMatrixView(this.ll, TriMatrixType.lower, false);
    }
    LT(): TriMatrixView {
        return new TriMatrixView(this.ll, TriMatrixType.upper, false);
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