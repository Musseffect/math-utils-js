import Matrix from "./denseMatrix";
import SparseMatrix from "./sparseMatrix";
import Triplet from "./triplet";
import { assert, swap } from "./utils";
import Vector from "./vector";

// row permutation - pre multiplied, col permutation - post multiplied
export default class PermutationMatrix {
    private permutations: number[];
    private _isRow: boolean;
    constructor(permutations: number[], isRow: boolean) {
        this.permutations = permutations.slice();
        this._isRow = isRow;
    }
    static identity(size: number, isRow: boolean): PermutationMatrix {
        let data = Array(size);
        for (let i = 0; i < size; ++i)
            data[i] = i;
        return new PermutationMatrix(data, isRow);
    }
    clone(): PermutationMatrix {
        return new PermutationMatrix(this.permutations.slice(), this._isRow);
    }
    isRow() {
        return this._isRow;
    }
    swap(i: number, j: number) {
        swap(this.permutations, i, j);
    }
    size() {
        return this.permutations.length;
    }
    findIndexByValue(value: number) {
        assert(value >= 0 && value < this.permutations.length, "Invalid value");
        let index = this.permutations.findIndex((v) => {
            return v == value;
        });
        assert(index != -1, "Invalid permutation");
        return index;
    }
    determinant(): number {
        let s = 0;
        for (let i = 0; i < this.permutations.length; ++i)
            s += Number(i == this.permutations[i]);
        return s & 1 ? -1 : 1;
    }
    isValid(): boolean {
        let values = new Array(this.permutations.length);
        values.fill(0);
        this.permutations.forEach((v) => {
            if (v >= 0 && v < this.permutations.length)
                values[v]++;
        });
        for (const value of values) {
            if (value != 1) return false;
        }
        return true;
    }
    static inverse(permutations: number[]): number[] {
        let values = new Array(permutations.length);
        permutations.forEach((v, i) => {
            values[v] = i;
        });
        return values;
    }
    value(i: number): number {
        return this.permutations[i];
    }
    inverse(): PermutationMatrix {
        return new PermutationMatrix(PermutationMatrix.inverse(this.permutations), this._isRow);
    }
    array(): number[] {
        return this.permutations;
    }
    permuteVector(v: Vector): Vector {
        let result = v.clone();
        for (let i = 0; i < result.size(); ++i)
            result.set(i, v.get(this.value(i)));
        return result;
    }
    permuteMatrix(m: Matrix): Matrix {
        let result = m.clone();
        for (let i = 0; i < result.numCols(); ++i) {
            for (let j = 0; j < result.numRows(); ++j) {
                let { row, column } = this.unpermuteIndex(j, i);
                result.set(j, i, m.get(row, column));
            }
        }
        return result;
    }
    permuteIndex(row: number, column: number) {
        if (this.isRow())
            row = this.findIndexByValue(row);
        else column = this.findIndexByValue(column);
        return { row, column };
    }
    unpermuteIndex(row: number, column: number) {
        if (this.isRow())
            row = this.permutations[row];
        else column = this.permutations[column];
        return { row, column };
    }
    get(row: number, column: number): number {
        if (this.isRow()) {
            if (column == this.permutations[row]) return 1;
        } else {
            if (row == this.permutations[column]) return 1;
        }
        return 0;
    }
    toMatrix(): Matrix {
        let m: Matrix = Matrix.empty(this.permutations.length, this.permutations.length);
        for (let i = 0; i < this.permutations.length; ++i) {
            const index = this.permutations[i];
            if (this.isRow())
                m.set(i, index, 1);
            else
                m.set(index, i, 1)
        }
        return m;
    }
    toTriplets(): Triplet[] {
        let triplets: Triplet[] = [];
        for (let i = 0; i < this.permutations.length; ++i) {
            const index = this.permutations[i];
            if (this.isRow())
                triplets.push({ row: i, column: index, value: 1 });
            else
                triplets.push({ row: index, column: i, value: 1 });
        }
        return triplets;
    }
    toSparseMatrix(): SparseMatrix {
        return SparseMatrix.fromTriplets(this.toTriplets(), this.permutations.length, this.permutations.length);
    }
}