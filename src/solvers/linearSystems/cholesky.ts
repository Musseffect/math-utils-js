//todo:

import matrix from "../../denseMatrix";
import { assert, near, SmallEpsilon } from "../../utils";
import vector from "../../vector";
import { InsufficientRankException } from "./exceptions";

const SolverName = "'cholesky'";

export default class cholesky {
    // Solve inplace
    static solve(A: matrix, b: vector, tolerance: number = SmallEpsilon) {
        assert(A.w == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.w == A.h, "Non-square matrix");

        let rank = b.size();
        for (let row = 0; row < rank; ++row) {
            for (let column = 0; column <= row; ++column) {
                let value = A.get(row, column);
                for (let k = 0; k < column; ++k)
                    value -= A.get(row, k) * A.get(column, k);
                if (row == column) {
                    if (value < 0.0) throw new InsufficientRankException(SolverName);
                    value = Math.sqrt(value);
                }
                else {
                    let denum = A.get(column, column);
                    if (near(denum, 0.0, tolerance)) throw new InsufficientRankException(SolverName);
                    value = value / A.get(column, column);
                }
                A.set(row, column, value);
                A.set(column, row, value);
            }
        }
        let y = vector.empty(rank);
        for (let row = 0; row < rank; ++row) {
            let value = b.get(row);
            for (let column = 0; column < row; ++column) {
                value -= A.get(row, column) * y.get(column);
            }
            value /= A.get(row, row);
            y.set(row, value);
        }
        let x = vector.empty(rank);
        for (let row = rank - 1; row >= 0; --row) {
            let value = y.get(row);
            for (let column = row + 1; column < rank; ++column) {
                value -= A.get(row, column) * x.get(column);
            }
            value /= A.get(row, row);
            x.set(row, value);
        }
        return x;
    }
}