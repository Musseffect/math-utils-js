import Matrix from "../../../denseMatrix";
import PermutationMatrix from "../../../permutationMatrix";
import { TriMatrixType, TriMatrixView } from "../../../triMatrixView";
import { assert, SmallEpsilon, swap } from "../../../utils";
import Vector from "../../../vector";
import { InsufficientRankException } from "./exceptions";

const SolverName = "'fullPivLU'";

/* LU decomposition with row and column permutations*/
export default class FullPivLU {
    protected lu: Matrix;
    protected p: PermutationMatrix;
    protected q: PermutationMatrix;
    protected A: Matrix;
    protected _rank: number;
    private compute(): void {
        this._rank = 0;
        throw new Error("Not implemented");
    }
    constructor(A: Matrix) {
        this.A = A;
        this.compute();
    }
    public rank(): number {
        return this._rank;
    }
    public L(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.lower, true);
    }
    public U(): TriMatrixView {
        return new TriMatrixView(this.lu, TriMatrixType.upper, false);
    }
    public Q(): Matrix {
        return this.q.toMatrix();
    }
    public P(): Matrix {
        return this.p.toMatrix();
    }
    static solveMatrix(A: Matrix, B: Matrix, tolerance: number = SmallEpsilon) {
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
    static solve(A: Matrix, b: Vector, tolerance: number = SmallEpsilon): Vector {
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