import matrix from "../../../denseMatrix";
import { assert, SmallEpsilon, swap } from "../../../utils";
import vector from "../../../vector";
import { InsufficientRankException } from "./exceptions";

const SolverName = "'partialPivLU'";

/* LU decomposition with row permutations*/
export default class PartialPivLU {
    static solve(A: matrix, b: vector, tolerance: number = SmallEpsilon) {
        assert(A.w == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.w == A.h, "Non-square matrix");
        let rank = b.size();
        let x = vector.empty(rank);
        let permutations = new Array(rank);
        for (let row = 0; row < rank; row++) {
            permutations[row] = row;
        }
        for (let step = 0; step < rank; step++) {
            let maxPivotValueRowIdx = step;
            let maxPivotValue = A.get(permutations[maxPivotValueRowIdx], step);
            for (let row = step + 1; row < rank; row++) {
                let value = A.get(permutations[row], step);
                if (Math.abs(value) > Math.abs(maxPivotValue)) {
                    maxPivotValueRowIdx = row;
                    maxPivotValue = value;
                }
            }

            if (Math.abs(maxPivotValue) < tolerance)
                throw new InsufficientRankException(SolverName);

            if (maxPivotValueRowIdx != step)
                swap(permutations, step, maxPivotValueRowIdx);

            for (let column = step + 1; column < rank; column++)
                A.set(permutations[step], column, A.get(permutations[step], column) / maxPivotValue);
            b.set(permutations[step], b.get(permutations[step]) / maxPivotValue);
            A.set(permutations[step], step, 1);

            for (let row = step + 1; row < rank; row++) {
                for (let column = step + 1; column < rank; column++)
                    A.set(permutations[row], column, A.get(permutations[row], column) - A.get(permutations[row], step) * A.get(permutations[step], column));
                b.set(permutations[row], b.get(permutations[row]) - A.get(permutations[row], step) * b.get(permutations[step]));
                A.set(permutations[row], step, 0);
            }
        }

        x.set(rank - 1, b.get(permutations[rank - 1]));
        for (let row = rank - 2; row > -1; row--) {
            let k = 0.;
            for (let column = row + 1; column < rank; column++)
                k = k + A.get(permutations[row], column) * x.get(column);

            x.set(row, b.get(permutations[row]) - k);
        }
        return x;
    }
}