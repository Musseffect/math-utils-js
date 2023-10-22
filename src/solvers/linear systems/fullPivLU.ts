import Matrix from "../../denseMatrix";
import PermutationMatrix from "../../permutationMatrix";
import { DiagonalType, TriMatrixType, TriMatrixView } from "../../triMatrixView";
import { assert, SmallestTolerance, SmallTolerance, swap } from "../../utils";
import Vector from "../../vector";
import { InsufficientRankException } from "./exceptions";

const SolverName = "'fullPivLU'";

/** LU decomposition with row and column permutations in the form PAQ = LU
 * 
 * L - lower unitriangular matrix
 * 
 * U - upper triangular matrix
 */
export default class FullPivLU {
    protected lu: Matrix = null;
    protected p: PermutationMatrix;
    protected q: PermutationMatrix;
    protected A: Matrix;
    protected _tolerance: number = SmallestTolerance;
    protected _rank: number;
    constructor(A: Matrix | null = null, tolerance: number = SmallestTolerance) {
        this.tolerance = tolerance;
        this.decompose(A);
    }
    set tolerance(value: number) {
        this._tolerance = value;
    }
    public rank(): number {
        return this._rank;
    }
    get L(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.lower, DiagonalType.Unit);
    }
    get U(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.upper, DiagonalType.Zero);
    }
    get LU(): Matrix {
        return this.lu;
    }
    get Q(): Matrix {
        return this.q.toMatrix();
    }
    get P(): Matrix {
        return this.p.toMatrix();
    }
    public decompose(A: Matrix | null): void {
        this.A = A;
        this.lu = null;
        this.p = null;
        this.q = null;
        if (A == null)
            return;
        assert(A.isSquare(), "Non-square matrix");
        this.p = PermutationMatrix.identity(this.A._numRows, true);
        this.q = PermutationMatrix.identity(this.A._numCols, false);
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
            if (Math.abs(maxPivot) < this._tolerance) {
                this._rank = step;
                this.lu = lu;
                return;
            }

            this.p.swap(step, maxPivotRow);
            this.q.swap(step, maxPivotColumn);

            lu.swapRows(step, maxPivotRow);
            lu.swapColumns(step, maxPivotColumn);
            // console.log(`Step ${step}`);
            // console.log(`rowPermutations ${this.p.toString()}, maxPivotRow ${maxPivotRow}`);
            // console.log(`columnPermutations ${this.q.toString()}, maxPivotColumn ${maxPivotColumn}`);
            // console.log(`Initial LU ${lu.toString()}`)
            /*const rowMat = this.p.toMatrix();
            const colMat = this.q.toMatrix();
            console.log(`Initial permuted LU ${Matrix.mul(Matrix.mul(rowMat, LU), colMat).toString()}`)
            console.log(`Initial Rhs ${Rhs.toString()}`)
            console.log(`Initial permuted Rhs ${Matrix.mul(Matrix.mul(rowMat, Rhs), colMat).toString()}`);
            */
            let ratio = lu.get(step, step) / maxPivot;
            for (let row = step + 1; row < lu._numRows; row++) {
                for (let column = step + 1; column < lu._numCols; column++)
                    lu.set(row, column, lu.get(row, column) - ratio * lu.get(row, step));
                lu.set(step, maxPivotColumn, ratio);
            }
            // console.log(`Result LU ${LU.toString()}`)
            // console.log(`Result permuted LU ${Matrix.mul(Matrix.mul(rowMat, LU), colMat).toString()}`)
            // console.log(`Result Rhs ${Rhs.toString()}`)
            // console.log(`Result permuted Rhs ${Matrix.mul(Matrix.mul(rowMat, Rhs), colMat).toString()}`);
        }
        lu = this.lu;
    }
    solver(rhs: Matrix | Vector): Vector {
        assert(this.A != null, "Decomposition is not available");
        if (rhs instanceof Matrix) {
            assert(this.A._numRows == rhs._numRows, "Incompatible RHS");
            throw new Error("Not implemented");

        } else {
            assert(this.A._numRows == rhs.size(), "Incompatible RHS");
            throw new Error("Not implemented");
        }
    }
    static solveMatrix(A: Matrix, B: Matrix, tolerance: number = SmallTolerance) {
        assert(A.height() == B.height(), "Not determined system");
        assert(A.isSquare(), "Non-square matrix");

        let LU = A.clone();
        let Rhs = B.clone();

        let rank = Rhs.height();
        let rowPermutations = new Array(rank);
        let columnPermutations = new Array(rank);
        for (let i = 0; i < rank; i++) {
            rowPermutations[i] = i;
            columnPermutations[i] = i;
        }
        let rowIdx = (idx: number) => rowPermutations[idx];
        let colIdx = (idx: number) => columnPermutations[idx];
        for (let step = 0; step < rank; step++) {
            // todo: extract this to a local function
            let maxPivotRow = step;
            let maxPivotColumn = step;
            let maxPivot = LU.get(rowIdx(maxPivotRow), colIdx(maxPivotColumn));
            for (let row = step; row < rank; ++row) {
                for (let column = step; column < rank; ++column) {
                    let value = LU.get(rowIdx(row), colIdx(column));
                    if (Math.abs(value) > Math.abs(maxPivot)) {
                        maxPivotRow = row;
                        maxPivotColumn = column;
                        maxPivot = value;
                    }
                }
            }
            if (Math.abs(maxPivot) < tolerance)
                throw new InsufficientRankException(SolverName, step);

            swap(rowPermutations, step, maxPivotRow);
            swap(columnPermutations, step, maxPivotColumn);

            console.log(`Step ${step}`);
            console.log(`rowPermutations ${rowPermutations}, maxPivotRow ${maxPivotRow}`);
            console.log(`columnPermutations ${columnPermutations}, maxPivotColumn ${maxPivotColumn}`);
            let stepRowIdx = rowIdx(step);

            console.log(`Initial LU ${LU.toString()}`)
            const rowMat = new PermutationMatrix(rowPermutations, true).toMatrix();
            const colMat = new PermutationMatrix(columnPermutations, false).toMatrix();
            console.log(`Initial permuted LU ${Matrix.mul(Matrix.mul(rowMat, LU), colMat).toString()}`)
            console.log(`Initial Rhs ${Rhs.toString()}`)
            console.log(`Initial permuted Rhs ${Matrix.mul(Matrix.mul(rowMat, Rhs), colMat).toString()}`);
            for (let row = step + 1; row < rank; row++) {
                let curRow = rowIdx(row);
                let ratio = LU.get(curRow, colIdx(step)) / maxPivot;
                for (let column = step + 1; column < rank; column++) {
                    let curColIdx = colIdx(column);
                    LU.set(curRow, curColIdx, LU.get(curRow, curColIdx) - ratio * LU.get(stepRowIdx, curColIdx));
                }
                for (let bCol = 0; bCol < Rhs.width(); ++bCol) {
                    Rhs.set(curRow, bCol, Rhs.get(curRow, bCol) - ratio * Rhs.get(stepRowIdx, bCol));
                }
                LU.set(curRow, colIdx(step), 0);
            }
            console.log(`Result LU ${LU.toString()}`)
            console.log(`Result permuted LU ${Matrix.mul(Matrix.mul(rowMat, LU), colMat).toString()}`)
            console.log(`Result Rhs ${Rhs.toString()}`)
            console.log(`Result permuted Rhs ${Matrix.mul(Matrix.mul(rowMat, Rhs), colMat).toString()}`);
        }

        let x = Matrix.empty(rank, Rhs.width());
        for (let bCol = 0; bCol < Rhs.width(); ++bCol) {
            let bColIdx = rowIdx(bCol);
            x.set(colIdx(rank - 1), bColIdx, Rhs.get(rowIdx(rank - 1), bColIdx) / LU.get(rowIdx(rank - 1), colIdx(rank - 1)));
            for (let row = rank - 2; row > -1; row--) {
                let curRowIdx = rowIdx(row);
                let k = 0.;
                for (let column = row + 1; column < rank; column++)
                    k = k + LU.get(curRowIdx, colIdx(column)) * x.get(colIdx(column), bColIdx);

                x.set(colIdx(row), bColIdx, (Rhs.get(curRowIdx, bColIdx) - k) / LU.get(curRowIdx, colIdx(row)));
            }
        }
        console.log(`X ${x.toString()}`);
        return x;
    }
    static solve(A: Matrix, b: Vector, tolerance: number = SmallTolerance): Vector {
        assert(A.width() == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.width() == A.height(), "Non-square matrix");

        let LU = A.clone();
        let rhs = b.clone();

        let rank = rhs.size();
        let x = Vector.empty(rank);
        let rowPermutations = new Array(rank);
        let columnPermutations = new Array(rank);
        for (let i = 0; i < rank; i++) {
            rowPermutations[i] = i;
            columnPermutations[i] = i;
        }
        let rowIdx = (idx: number) => rowPermutations[idx];
        let colIdx = (idx: number) => columnPermutations[idx];
        for (let step = 0; step < rank; step++) {
            let maxPivotRow = step;
            let maxPivotColumn = step;
            let maxPivot = LU.get(rowIdx(maxPivotRow), colIdx(maxPivotColumn));
            for (let row = step; row < rank; ++row) {
                for (let column = step; column < rank; ++column) {
                    let value = LU.get(rowIdx(row), colIdx(column));
                    if (Math.abs(value) > Math.abs(maxPivot)) {
                        maxPivotRow = row;
                        maxPivotColumn = column;
                        maxPivot = value;
                    }
                }
            }
            if (Math.abs(maxPivot) < tolerance)
                throw new InsufficientRankException(SolverName, step);

            swap(rowPermutations, step, maxPivotRow);
            swap(columnPermutations, step, maxPivotColumn);

            for (let column = step + 1; column < rank; column++)
                LU.set(rowIdx(step), colIdx(column), LU.get(rowIdx(step), colIdx(column)) / maxPivot);
            rhs.set(rowIdx(step), rhs.get(rowIdx(step)) / maxPivot);
            LU.set(rowIdx(step), colIdx(step), 1);

            for (let row = step + 1; row < rank; row++) {
                let firstValue = LU.get(rowIdx(row), colIdx(step));
                for (let column = step + 1; column < rank; column++)
                    LU.set(rowIdx(row), colIdx(column), LU.get(rowIdx(row), colIdx(column)) - firstValue * LU.get(rowIdx(step), colIdx(column)));
                rhs.set(rowIdx(row), rhs.get(rowIdx(row)) - firstValue * rhs.get(rowIdx(step)));
                LU.set(rowIdx(row), colIdx(step), 0);
            }
        }

        x.set(colIdx(rank - 1), rhs.get(rowIdx(rank - 1)));
        for (let row = rank - 2; row > -1; row--) {
            let value = rhs.get(rowIdx(row));
            for (let column = row + 1; column < rank; column++)
                value -= LU.get(rowIdx(row), colIdx(column)) * x.get(colIdx(column));

            x.set(colIdx(row), value);
        }
        return x;
    }
}