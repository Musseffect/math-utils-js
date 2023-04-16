import Matrix from "../../../denseMatrix";
import PermutationMatrix from "../../../permutationMatrix";
import SparseMatrix from "../../../SparseMatrix";
import { assert, SmallEpsilon, swap } from "../../../utils";
import Vector from "../../../vector";
import { TriMatrixType, TriMatrixView } from "../../../triMatrixView";

const SolverName = "'PartialPivLU'";

/* 
    LU decomposition with row permutations in the form PA = LU
    for Ax = b square problems
*/
export default class PartialPivLU {
    private lu: Matrix;
    private p: PermutationMatrix;
    private A: Matrix;
    private compute(): void {
        throw new Error("Not implemented");
    }
    constructor(A: Matrix, tolerance: number = SmallEpsilon) {
        this.A = A;
        this.compute();
    }
    static solveMatrix(A: Matrix, B: Matrix, tolerance: number = SmallEpsilon): Matrix {
        assert(A.height() == B.height(), "Not determined system");
        assert(A.isSquare(), "Non-square matrix");
        let rank = B.height();

        let LU = A.clone();
        let Rhs = B.clone();

        let permutations = new Array(rank);
        for (let row = 0; row < rank; row++)
            permutations[row] = row;

        for (let step = 0; step < rank; step++) {
            let maxPivotRow = step;
            let maxPivotValue = LU.get(permutations[maxPivotRow], step);
            for (let row = step + 1; row < rank; row++) {
                let value = LU.get(permutations[row], step);
                if (Math.abs(value) > Math.abs(maxPivotValue)) {
                    maxPivotRow = row;
                    maxPivotValue = value;
                }
            }

            swap(permutations, step, maxPivotRow);

            let stepRow = permutations[step];

            for (let row = step + 1; row < rank; row++) {
                let curRow = permutations[row];
                let ratio = LU.get(curRow, step) / maxPivotValue;
                for (let column = step + 1; column < rank; column++)
                    LU.set(curRow, column, LU.get(curRow, column) - ratio * LU.get(stepRow, column));
                for (let bCol = 0; bCol < Rhs.width(); ++bCol)
                    Rhs.set(curRow, bCol, Rhs.get(curRow, bCol) - ratio * Rhs.get(stepRow, bCol));
                LU.set(curRow, step, 0);
            }
        }
        let x = Matrix.empty(rank, Rhs.width());
        for (let bCol = 0; bCol < Rhs.width(); ++bCol) {
            x.set(rank - 1, bCol, Rhs.get(permutations[rank - 1], bCol) / LU.get(permutations[rank - 1], rank - 1));
            for (let row = rank - 2; row > -1; row--) {
                let curRow = permutations[row];
                let k = 0.;
                for (let column = row + 1; column < rank; column++)
                    k = k + LU.get(curRow, column) * x.get(column, bCol);

                x.set(row, bCol, (Rhs.get(curRow, bCol) - k) / LU.get(curRow, row));
            }
        }
        return x;
    }
    static solve(A: Matrix, b: Vector, tolerance: number = SmallEpsilon): Vector {
        assert(A.width() == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.width() == A.height(), "Non-square matrix");

        let LU = A.clone();
        let rhs = b.clone();

        let rank = rhs.size();
        let permutations = new Array(rank);
        for (let row = 0; row < rank; row++)
            permutations[row] = row;

        for (let step = 0; step < rank; step++) {
            let maxPivotRow = step;
            let maxPivotValue = LU.get(permutations[maxPivotRow], step);
            for (let row = step + 1; row < rank; row++) {
                let value = LU.get(permutations[row], step);
                if (Math.abs(value) > Math.abs(maxPivotValue)) {
                    maxPivotRow = row;
                    maxPivotValue = value;
                }
            }

            swap(permutations, step, maxPivotRow);

            let stepRow = permutations[step];

            for (let row = step + 1; row < rank; row++) {
                let curRow = permutations[row];
                let ratio = LU.get(curRow, step) / maxPivotValue;
                for (let column = step + 1; column < rank; column++)
                    LU.set(curRow, column, LU.get(curRow, column) - ratio * LU.get(stepRow, column));
                rhs.set(curRow, rhs.get(curRow) - ratio * rhs.get(stepRow));
                LU.set(curRow, step, 0);
            }
        }
        let x = Vector.empty(rank);
        x.set(rank - 1, rhs.get(permutations[rank - 1]) / LU.get(permutations[rank - 1], rank - 1));
        for (let row = rank - 2; row > -1; row--) {
            let curRow = permutations[row];
            let k = 0.;
            for (let column = row + 1; column < rank; column++)
                k = k + LU.get(curRow, column) * x.get(column);

            x.set(row, (rhs.get(curRow) - k) / LU.get(curRow, row));
        }
        return x;
    }
    public L(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.lower, true);
    }
    public U(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.upper, false);
    }
    public P(): Matrix {
        return this.p.toMatrix();
    }
}