import Matrix from "../../denseMatrix";
import PermutationMatrix from "../../permutationMatrix";
import { assert, SmallestTolerance, SmallTolerance, swap } from "../../utils";
import Vector from "../../vector";
import { DiagonalType, TriMatrixType, TriMatrixView } from "../../triMatrixView";

const SolverName = "'PartialPivLU'";

/* 
    LU decomposition with row permutations in the form PA = LU
    for Ax = b square problems
*/
export default class PartialPivLU {
    private lu: Matrix = null;
    private p: PermutationMatrix;
    private A: Matrix;
    private _tolerance: number = SmallestTolerance;
    private decompose(A: Matrix | null): void {
        this.A = A;
        this.lu = null;
        this.p = null;
        if (A == null)
            return;
        assert(A.isSquare(), "Non-square matrix");
        this.p = PermutationMatrix.identity(this.A._numRows, true);
        let lu: Matrix = this.A.clone();
        // todo: check for rectangular matrices
        for (let step = 0; step < this.lu._numRows; step++) {
            let maxPivotRow = step;
            let maxPivotColumn = step;
            let maxPivot = lu.get(step, step);
            for (let row = step; row < lu._numRows; ++row) {
                for (let column = step; column < lu._numCols; ++column) {
                    let value = lu.get(row, column);
                    if (Math.abs(value) > Math.abs(maxPivot)) {
                        maxPivotRow = row;
                        maxPivotColumn = column;
                        maxPivot = value;
                    }
                }
            }
            if (Math.abs(maxPivot) < this._tolerance)
                return;

            this.p.swap(step, maxPivotRow);

            lu.swapRows(step, maxPivotRow);
            let ratio = lu.get(step, step) / maxPivot;
            for (let row = step + 1; row < lu._numRows; row++) {
                for (let column = step + 1; column < lu._numCols; column++)
                    lu.set(row, column, lu.get(row, column) - ratio * lu.get(row, step));
                lu.set(step, maxPivotColumn, ratio);
            }
        }
        this.lu = lu;
    }
    get L(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.lower, DiagonalType.Unit);
    }
    get U(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.upper, DiagonalType.Existing);
    }
    get LU(): Matrix {
        return this.lu;
    }
    get P(): PermutationMatrix {
        return this.p;
    }
    set tolerance(value: number) {
        this._tolerance = value;
    }
    constructor(A: Matrix | null = null, tolerance: number = SmallTolerance) {
        this.tolerance = tolerance;
        this.A = A;
        this.decompose(A);
    }
    static solveMatrix(A: Matrix, B: Matrix, tolerance: number = SmallTolerance): Matrix {
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
    isSingular(): boolean { return this.lu == null; }
    determinant(): number {
        if (this.isSingular()) return 0.0;
        let determinant = this.p.determinant();
        for (let i = 0; i < this.lu.width(); ++i)
            determinant *= this.lu.get(i, i);
        return determinant;
    }
    inverse(): Matrix | null {
        if (this.lu == null) return null;
        throw new Error("Not implemented");
    }
    static solve(A: Matrix, b: Vector, tolerance: number = SmallTolerance): Vector {
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
}