import matrix from "../../denseMatrix";
import { assert, SmallEpsilon, swap } from "../../utils";
import vector from "../../vector";
import { InsufficientRankException } from "./exceptions";


/* LU decomposition with row and column permutations*/
export default class FullPivLU {
    static solve(A: matrix, b: vector, tolerance: number = SmallEpsilon) {
        assert(A.w == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.w == A.h, "Non-square matrix");

        let rank = b.size();
        let x = vector.empty(rank);
        let rowPermutations = new Array(rank);
        let columnPermutations = new Array(rank);
        for (let i = 0; i < rank; i++) {
            rowPermutations[i] = i;
            columnPermutations[i] = i;
        }
        let rowIdx = (idx: number) => rowPermutations[idx];
        let colIdx = (idx: number) => columnPermutations[idx];
        for (let step = 0; step < rank; step++) {
            let maxPivotRowIdx = step;
            let maxPivotColumnIdx = step;
            let maxPivot = A.get(rowIdx(maxPivotRowIdx), colIdx(maxPivotColumnIdx));
            for (let row = step; row < rank; ++row) {
                for (let column = step; column < rank; ++column) {
                    let value = A.get(rowIdx(row), colIdx(column));
                    if (Math.abs(value) > Math.abs(maxPivot)) {
                        maxPivotRowIdx = row;
                        maxPivotColumnIdx = column;
                        maxPivot = value;
                    }
                }
            }

            if (Math.abs(maxPivot) < tolerance)
                throw new InsufficientRankException("'fullPivLU'");

            if (maxPivotRowIdx != step)
                swap(rowPermutations, step, maxPivotRowIdx);
            if (maxPivotColumnIdx != step)
                swap(columnPermutations, step, maxPivotColumnIdx);

            for (let column = step + 1; column < rank; column++)
                A.set(rowIdx(step), colIdx(column), A.get(rowIdx(step), colIdx(column)) / maxPivot);
            b.set(rowIdx(step), b.get(rowIdx(step)) / maxPivot);
            A.set(rowIdx(step), colIdx(step), 1);

            for (let row = step + 1; row < rank; row++) {
                let firstValue = A.get(rowIdx(row), colIdx(step));
                for (let column = step + 1; column < rank; column++)
                    A.set(rowIdx(row), colIdx(column), A.get(rowIdx(row), colIdx(column)) - firstValue * A.get(rowIdx(step), colIdx(column)));
                b.set(rowIdx(row), b.get(rowIdx(row)) - firstValue * b.get(rowIdx(step)));
                A.set(rowIdx(row), colIdx(step), 0);
            }
        }

        x.set(colIdx(rank - 1), b.get(rowIdx(rank - 1)));
        for (let row = rank - 2; row > -1; row--) {
            let value = b.get(rowIdx(row));
            for (let column = row + 1; column < rank; column++)
                value -= A.get(rowIdx(row), colIdx(column)) * x.get(colIdx(column));

            x.set(colIdx(row), value);
        }
        return x;
    }
}