import { assert, Tolerance } from "./utils";


export default abstract class AbstractDenseMatrix {
    data: number[];
    constructor(data: number[]) {
        this.data = data;
    }
    abstract numCols(): number;
    abstract numRows(): number;
    get(row: number, col: number): number {
        return this.data[row * this.numCols() + col];
    }
    set(row: number, col: number, value: number): void {
        this.data[row * this.numCols() + col] = value;
    }
    values(): number[] {
        return this.data;
    }
    addSelf(a: AbstractDenseMatrix): void {
        const rows = this.numRows();
        const cols = this.numCols();
        assert(a.numRows() == rows && a.numCols() == cols, "Sizes don't match");
        for (let i = 0; i < this.data.length; ++i)
            this.data[i] += a.data[i];
    }
    subSelf(a: AbstractDenseMatrix): void {
        const rows = this.numRows();
        const cols = this.numCols();
        assert(a.numRows() == rows && a.numCols() == cols, "Sizes don't match");
        for (let i = 0; i < this.data.length; ++i)
            this.data[i] -= a.data[i];
    }
    isSquare(): boolean {
        return this.numCols() == this.numRows();
    }
    mulSelfPointwise(a: AbstractDenseMatrix): void {
        const rows = this.numRows();
        const cols = this.numCols();
        assert(a.numRows() == rows && a.numCols() == cols, "Sizes don't match");
        for (let i = 0; i < this.data.length; ++i)
            this.data[i] *= a.data[i];
    }
    divSelfPointwise(a: AbstractDenseMatrix): void {
        const rows = this.numRows();
        const cols = this.numCols();
        assert(a.numRows() == rows && a.numCols() == cols, "Sizes don't match");
        for (let i = 0; i < this.data.length; ++i)
            this.data[i] *= a.data[i];
    }
    scaleSelf(scalar: number): void {
        for (let i = 0; i < this.data.length; ++i)
            this.data[i] *= scalar;
    }
    l2SquaredNorm(): number {
        let result = 0.0;
        for (let value of this.data)
            result += value * value;
        return result;
    }
    l2Norm(): number {
        return Math.sqrt(this.l2SquaredNorm());
    }
    lInfNorm(): number {
        let result = 0.0;
        for (let value of this.data)
            result = Math.max(result, Math.abs(value));
        return result;
    }
}