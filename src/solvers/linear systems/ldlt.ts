import Matrix from "../../denseMatrix";
import { TriMatrixView } from "../../triMatrixView";
import { assert } from "../../utils";
import Vector from "../../vector";


class TriangularMatrix {
    isLower: boolean;
    data: number[];
    size: number;
    constructor(size: number, isLower: boolean = true) {
        this.size = size;
        this.data = new Array(size * (size + 1) / 2).fill(0);
        this.isLower = isLower;
    }
    get(row: number, col: number) {
        let i = this.isLower ? row : col;
        let j = this.isLower ? col : row;
        if (i >= j) return 0.0;
        return this.data[(i * (i + 1)) / 2 + j];
    }
    set(row: number, col: number, value: number) {
        let i = this.isLower ? row : col;
        let j = this.isLower ? col : row;
        assert(i >= j, "Invalid index");
        this.data[(i * (i + 1)) / 2 + j] = value;
    }
}

export default class LDLT {
    LDLT: Matrix;
    constructor(tolerance: number) {
    }
    factorize(A: Matrix) {
        throw new Error("Not implemented");
    }
    solve(rhs: Vector): Vector {
        throw new Error("Not implemented");
    }
    L(): TriMatrixView {
        throw new Error("Not implemented");
    }
    D(): Vector {
        throw new Error("Not implemented");
    }
    inverse(): Matrix {
        assert(this.LDLT != null, "Factorization is not available");
        throw new Error("Not implemented");
    }
    static solve(A: Matrix, rhs: Vector) {
        const rank = A.width();
        let L = new TriangularMatrix(rank, true);
        for (let row = 0; row < rank; ++row) {
            for (let col = 0; col < row; ++col) {
                let value = A.get(row, col);
                for (let i = 0; i < col; ++i) {
                    value -= L.get(row, i) * L.get(col, i) * L.get(i, i);
                }
                L.set(row, col, value / L.get(col, col));
            }
            let value = A.get(row, row);
            for (let i = 0; i < row; ++i)
                value -= L.get(row, i) * L.get(row, i) * L.get(row, row);
            L.set(row, row, value);
        }
        let y = Vector.empty(rank);
        for (let row = 0; row < rank; ++row) {
            let value = rhs.get(row);
            for (let column = 0; column < row; ++column)
                value -= L.get(row, column) * y.get(column);
            y.set(row, value / L.get(row, row));
        }
        let x = Vector.empty(rank);
        for (let row = rank - 1; row >= 0; --row) {
            let value = y.get(row);
            for (let column = row + 1; column < rank; ++column)
                value -= L.get(column, row) * x.get(column);
            x.set(row, value);
        }
    }
}