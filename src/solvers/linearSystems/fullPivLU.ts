import matrix from "../../denseMatrix";
import { assert, SmallEpsilon, swap } from "../../utils";
import vector from "../../vector";


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
        for (let step = 0; step < rank; step++) {
            let maxPivotValueRowIdx = step;
            let maxPivotValueColumnIdx = step;
            let maxPivotValue = A.get(rowPermutations[maxPivotValueRowIdx], columnPermutations[maxPivotValueColumnIdx]);
            for (let row = step; row < rank; ++row) {
                for (let column = step; column < rank; ++column) {
                    let value = A.get(rowPermutations[row], columnPermutations[column]);
                    if (Math.abs(value) > Math.abs(maxPivotValue)) {
                        maxPivotValueRowIdx = row;
                        maxPivotValueColumnIdx = column;
                        maxPivotValue = value;
                    }
                }
            }

            if (Math.abs(maxPivotValue) < tolerance)
                throw new InsufficientRankException();

            if (maxPivotValueRowIdx != step)
                swap(rowPermutations, step, maxPivotValueRowIdx);
            if (maxPivotValueColumnIdx != step)
                swap(columnPermutations, step, maxPivotValueColumnIdx);

            for (let column = step + 1; column < rank; column++)
                A.set(rowPermutations[step], columnPermutations[column], A.get(rowPermutations[step], columnPermutations[column]) / maxPivotValue);
            b.set(rowPermutations[step], b.get(rowPermutations[step]) / maxPivotValue);
            A.set(rowPermutations[step], columnPermutations[step], 1);

            for (let row = step + 1; row < rank; row++) {
                for (let column = step + 1; column < rank; column++)
                    A.set(rowPermutations[row], columnPermutations[column], A.get(rowPermutations[row], columnPermutations[column]) - A.get(rowPermutations[row], columnPermutations[step]) * A.get(rowPermutations[step], columnPermutations[column]));
                b.set(rowPermutations[row], b.get(rowPermutations[row]) - A.get(rowPermutations[row], columnPermutations[step]) * b.get(columnPermutations[step]));
                A.set(rowPermutations[row], columnPermutations[step], 0);
            }
        }

        x.set(rank - 1, b.get(rowPermutations[rank - 1]));
        for (let row = rank - 2; row > -1; row--) {
            let k = 0.;
            for (let column = row + 1; column < rank; column++)
                k = k + A.get(rowPermutations[row], columnPermutations[column]) * x.get(columnPermutations[column]);

            x.set(row, b.get(rowPermutations[row]) - k);
        }
        return x;
    }
}