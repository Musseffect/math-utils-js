import matrix from "./denseMatrix";
import sparseVector from "./sparseVector";
import { assert, SmallestEpsilon } from "./utils";
import vector from "./vector";

export interface triplet {
    row: number,
    column: number,
    value: number
}

export class CellRef {
    private state: number;
    private matrix: sparseMatrix;
    private idx: number;
    constructor(matrix: sparseMatrix, idx: number, state: number) {
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
export default class sparseMatrix {
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
    clone(): sparseMatrix {
        let result = new sparseMatrix(this.numRows, this.numCols, this.tolerance);
        result.innerIndices = this.innerIndices.slice();
        result.nonZeroElements = this.nonZeroElements.slice();
        result.outerStarts = this.outerStarts.slice();
        return result;
    }
    static identity(): sparseMatrix {
        throw new Error("Not implemented");
    }
    static add(m1: sparseMatrix, m2: sparseMatrix): sparseMatrix {
        throw new Error("Not implemented");
    }
    static sub(m1: sparseMatrix, m2: sparseMatrix): sparseMatrix {
        throw new Error("Not implemented");
    }
    static mul(m1: sparseMatrix, m2: sparseMatrix): sparseMatrix {
        throw new Error("Not implemented");
    }
    determinant(): number {
        throw new Error("Not implemented");
    }
    static postMulSparse(m: sparseMatrix, v: sparseVector): sparseVector {
        throw new Error("Not implemented");
    }
    static preMulSparse(v: sparseVector, m: sparseMatrix): sparseVector {
        throw new Error("Not implemented");
    }
    static postMul(m: sparseMatrix, v: vector): vector {
        throw new Error("Not implemented");
    }
    static preMul(v: vector, m: sparseMatrix): vector {
        throw new Error("Not implemented");
    }
    scale(scalar: number): sparseMatrix {
        return this;
    }
    static scale(m: sparseMatrix, scalar: number): sparseMatrix {
        return m.clone().scale(scalar);
    }
    rowVector(row: number): sparseVector {
        throw new Error("Not implemented");
    }
    columnVector(column: number): sparseVector {
        throw new Error("Not implemented");
    }
    protected advanceState(): void {
        ++this.state;
    }
    reserve(rows: number, cols: number) {
        throw new Error("Not implemented");
    }
    static fromDense(dense: matrix): sparseMatrix {
        throw new Error("Not implemented");
    }
    static fromTriplets(triplets: triplet[], numRows: number, numCols: number, tolerance: number = SmallestEpsilon): sparseMatrix {
        // sorted in ascending "row by row" order
        triplets.sort((a: triplet, b: triplet) => {
            let rowSign = a.row - b.row;
            if (rowSign != 0) return rowSign;
            return a.column - b.column;
        });
        let result = new sparseMatrix(numRows, numCols, tolerance);
        if (triplets.length == 0) return result;
        let currentRow = triplets[0].row;
        for (let i = 0; i < triplets.length; ++i) {
            assert(triplets[i].row <= numRows && triplets[i].row >= 0, "Invalid row index");
            assert(triplets[i].column <= numCols && triplets[i].column >= 0, "Invalid column index");
            if (triplets[i].row != currentRow) {
                for (let row = currentRow; row < triplets[i].row; ++row)
                    result.outerStarts.push(result.nonZeroElements.length);
                currentRow = triplets[i].row;
            }
            result.nonZeroElements.push(triplets[i].value);
            result.innerIndices.push(triplets[i].column);
        }
        throw new Error("Not implemented");
    }
    compress() {
        this.advanceState();
        throw new Error("Not implemented");
    }
    height(): number { return this.numRows; }
    width(): number { return this.numCols; }
    numNonZeroes(): number {
        throw new Error("Not implemented");
    }
    coeffRef(row: number, column: number): CellRef {
        //find value if exist
        // create value if not
        // return ref
        throw new Error("Not implemented");
    }
    insert(row: number, column: number, value: number) {
        throw new Error("Not implemented");
    }
    set(row: number, column: number, value: number) {
        throw new Error("Not implemented");
    }
    get(row: number, column: number): number {
        throw new Error("Not implemented");
    }
}