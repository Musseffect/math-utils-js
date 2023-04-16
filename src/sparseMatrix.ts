import Matrix from "./denseMatrix";
import sparseVector from "./sparseVector";
import Triplet from "./triplet";
import { assert, SmallEpsilon, SmallestEpsilon } from "./utils";
import vector from "./vector";

export class CellRef {
    private state: number;
    private matrix: SparseMatrix;
    private idx: number;
    constructor(matrix: SparseMatrix, idx: number, state: number) {
        this.matrix = matrix;
        this.idx = idx;
        this.state = state;
    }
    get(): number {
        assert(this.state == this.matrix.state, "Invalid ref");
        return this.matrix.nonZeroElements[this.idx];
    }
    set(value: number) {
        assert(this.state == this.matrix.state, "Invalid ref");
        this.matrix.nonZeroElements[this.idx] = value;
    }
}
// CSR
export default class SparseMatrix {
    state: number = 0;
    // column indices sizeof(NNZ)
    innerIndices: number[];
    // values sizeof(NNZ)
    nonZeroElements: number[];
    // row starts
    outerStarts: number[];

    numRows: number;
    numCols: number;
    tolerance: number;
    constructor(numRows: number, numCols: number, tolerance: number = SmallestEpsilon) {
        this.numRows = numRows;
        this.numCols = numCols;
        this.innerIndices = [];
        this.nonZeroElements = [];
        this.outerStarts = [];
    }
    clone(): SparseMatrix {
        let result = new SparseMatrix(this.numRows, this.numCols, this.tolerance);
        result.innerIndices = this.innerIndices.slice();
        result.nonZeroElements = this.nonZeroElements.slice();
        result.outerStarts = this.outerStarts.slice();
        return result;
    }
    l2Norm(): number {
        let result = 0.0;
        for (let value of this.nonZeroElements)
            result += value * value;
        return Math.sqrt(result);
    }
    lInfNorm(): number {
        let result = 0.0;
        for (let value of this.nonZeroElements)
            result = Math.max(result, Math.abs(value));
        return result;
    }
    static near(m1: SparseMatrix, m2: SparseMatrix, tolerance: number = SmallEpsilon): boolean {
        assert(m1.width() == m2.width() && m1.height() == m2.height(), "Incompatible sizes");
        throw new Error("Method not implemented.");
    }
    static identity(size: number): SparseMatrix {
        throw new Error("Method not implemented.");
    }
    static add(m1: SparseMatrix, m2: SparseMatrix): SparseMatrix {
        throw new Error("Method not implemented.");
    }
    static sub(m1: SparseMatrix, m2: SparseMatrix): SparseMatrix {
        throw new Error("Method not implemented.");
    }
    static mul(m1: SparseMatrix, m2: SparseMatrix): SparseMatrix {
        throw new Error("Method not implemented.");
    }
    determinant(): number {
        throw new Error("Method not implemented.");
    }
    static postMulSparse(m: SparseMatrix, v: sparseVector): sparseVector {
        throw new Error("Method not implemented.");
    }
    static preMulSparse(v: sparseVector, m: SparseMatrix): sparseVector {
        throw new Error("Method not implemented.");
    }
    static postMul(m: SparseMatrix, v: vector): vector {
        throw new Error("Method not implemented.");
    }
    static preMul(v: vector, m: SparseMatrix): vector {
        throw new Error("Method not implemented.");
    }
    scale(scalar: number): SparseMatrix {
        for (let i = 0; i < this.nonZeroElements.length; ++i)
            this.nonZeroElements[i] *= scalar;
        return this;
    }
    static scale(m: SparseMatrix, scalar: number): SparseMatrix {
        return m.clone().scale(scalar);
    }
    rowVector(row: number): sparseVector {
        throw new Error("Method not implemented.");
    }
    columnVector(column: number): sparseVector {
        throw new Error("Method not implemented.");
    }
    protected advanceState(): void {
        ++this.state;
    }
    reserve(rows: number, cols: number) {
        throw new Error("Method not implemented.");
    }
    static fromDense(dense: Matrix, tolerance: number): SparseMatrix {
        throw new Error("Method not implemented.");
    }
    static fromTriplets(triplets: Triplet[], numRows: number, numCols: number, tolerance: number = SmallestEpsilon): SparseMatrix {
        // sorted in ascending "row by row" order
        triplets.sort((a: Triplet, b: Triplet) => {
            let rowSign = a.row - b.row;
            if (rowSign != 0) return rowSign;
            return a.column - b.column;
        });
        let result = new SparseMatrix(numRows, numCols, tolerance);
        if (triplets.length == 0) return result;
        result.outerStarts.push(0);
        let currentRow = triplets[0].row;
        for (let i = 0; i < triplets.length; ++i) {
            assert(triplets[i].row <= numRows && triplets[i].row >= 0, "Invalid row index");
            assert(triplets[i].column <= numCols && triplets[i].column >= 0, "Invalid column index");
            if (triplets[i].row != currentRow) {
                for (let row = currentRow + 1; row <= triplets[i].row; ++row)
                    result.outerStarts.push(result.nonZeroElements.length);
                currentRow = triplets[i].row;
            } else if (i > 0 && triplets[i].column == triplets[i - 1].column) {
                result.nonZeroElements[result.nonZeroElements.length - 1] = triplets[i].value;
            }
            result.nonZeroElements.push(triplets[i].value);
            result.innerIndices.push(triplets[i].column);
        }
        for (let row = currentRow + 1; row <= numRows; ++row)
            result.outerStarts.push(result.nonZeroElements.length);
        assert(result.outerStarts.length == numRows + 1, "result.outerStarts.length == numRows + 1");
        //throw new Error("Method not implemented.");
        return result;
    }
    toDense(): Matrix {
        let result = Matrix.empty(this.numRows, this.numCols);
        for (let rowIdx = 0; rowIdx < this.numRows; ++rowIdx) {
            let startIdx = this.outerStarts[rowIdx];
            let endIdx = this.outerStarts[rowIdx + 1];
            for (let j = startIdx; j < endIdx; ++j) {
                let colIdx = this.innerIndices[j];
                let value = this.nonZeroElements[j];
                result.set(rowIdx, colIdx, value);
            }
        }
        return result;
    }
    compress() {
        this.advanceState();
        throw new Error("Method not implemented.");
    }
    height(): number { return this.numRows; }
    width(): number { return this.numCols; }
    numNonZeroes(): number {
        throw new Error("Method not implemented.");
    }
    coeffRef(row: number, column: number): CellRef {
        //find value if exist
        // create value if not
        // return ref
        throw new Error("Method not implemented.");
    }
    insert(row: number, column: number, value: number) {
        throw new Error("Method not implemented.");
    }
    set(row: number, column: number, value: number) {
        throw new Error("Method not implemented.");
    }
    get(row: number, column: number): number {
        throw new Error("Method not implemented.");
    }
    inverse(): any {
        throw new Error("Method not implemented.");
    }
}