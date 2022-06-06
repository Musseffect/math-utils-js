import sparseVector from "./sparseVector";
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
    tolerance: number;
    constructor(rows: number, cols: number, rowMajor: boolean = true, tolerance: number) {
        this.rows = rows;
        this.cols = cols;
        this.rowMajor = rowMajor;
    }
    clone(): sparseMatrix {

    }
    static identity(): sparseMatrix {

    }
    static add(m1: sparseMatrix, m2: sparseMatrix): sparseMatrix {

    }
    static sub(m1: sparseMatrix, m2: sparseMatrix): sparseMatrix {

    }
    static mul(m1: sparseMatrix, m2: sparseMatrix): sparseMatrix {

    }
    determinant(): number {

    }
    static postMul(m: sparseMatrix, v: sparseVector): sparseVector {

    }
    static preMul(v: sparseVector, m: sparseMatrix): sparseVector {

    }
    scale(scalar: number): sparseMatrix {
        return this;
    }
    static scale(scalar: number): sparseMatrix {
        return this.clone().scale(scalar);
    }
    rowVector(row: number): sparseVector {

    }
    columnVector(column: number): sparseVector {

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
    get(row: number, column: number): number {
        throw new Error("Not implemented");
    }
}