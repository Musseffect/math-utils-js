import Matrix from "./denseMatrix";
import SparseMatrix from "./sparseMatrix";
import Triplet from "./triplet";
import { swap } from "./utils";

// row permutation - pre multiplied, col permutation - post multiplied
export default class PermutationMatrix {
    permutations: number[];
    isRow: boolean;
    constructor(permutations: number[], isRow: boolean) {
        this.permutations = permutations.slice();
        this.isRow = isRow;
    }
    swap(i: number, j: number) {
        swap(this.permutations, i, j);
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
            values[i] = v;
        });
        return values;
    }
    value(i: number): number {
        return this.permutations[i];
    }
    inverse(): PermutationMatrix {
        return new PermutationMatrix(PermutationMatrix.inverse(this.permutations), this.isRow);
    }
    permutIndex(row: number, column: number) {
        if (this.isRow)
            row = this.permutations[row];
        else column = this.permutations[column];
        return { row, column };
    }
    get(row: number, column: number): number {
        if (this.isRow) {
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
            if (this.isRow)
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
            if (this.isRow)
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