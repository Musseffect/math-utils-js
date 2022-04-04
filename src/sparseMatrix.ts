import { assert } from "./utils";

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
    }
    get(): number {
        assert(this.state == this.matrix.state, "Invalid ref");
        return this.matrix.values[this.idx];
    }
    set(value: number) {
        assert(this.state == this.matrix.state, "Invalid ref");
        this.matrix.values[this.idx] = value;
    }
}

export default class sparseMatrix {
    state: number = 0;
    rowMajor: boolean = true;
    values: number[];
    innerIndices: number[];
    outerStarts: number[];
    innerNonZeroes: number[];
    rows: number;
    cols: number;
    constructor(rows: number, cols: number, rowMajor: boolean = true) {
        this.rows = rows;
        this.cols = cols;
        this.rowMajor = rowMajor;
    }
    protected advanceState(): void {
        ++this.state;
    }
    reserve(rows: number, cols: number) {
        throw new Error("Not implemented");
    }
    fromTriplets(triplets: triplet[]): sparseMatrix {
        throw new Error("Not implemented");
    }
    compress() {
        this.advanceState();
        throw new Error("Not implemented");
    }
    numRows(): number { return this.rows; }
    numCols(): number { return this.cols; }
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
    read(row: number, column: number): number {
        throw new Error("Not implemented");
    }
}