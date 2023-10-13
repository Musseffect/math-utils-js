import mat from "./abstractDenseMatrix";
import FullPivLU from "./solvers/linear systems/fullPivLU";
import Triplet from "./triplet";
import { assert, near, SmallTolerance, SmallestTolerance, swap } from "./utils";
import vector from "./vector";

export default class Matrix extends mat {
    _numCols: number;
    _numRows: number;
    data: number[];
    constructor(data: number[], numRows: number, numCols: number) {
        super(data);
        this.data = data;
        this._numCols = numCols;
        this._numRows = numRows;
        assert(this.data.length == numRows * numCols, "Wrong size of data array");
    }
    numCols(): number {
        return this._numCols;
    }
    numRows(): number {
        return this._numRows;
    }
    isSymmetric(tolerance: number = SmallTolerance): boolean {
        if (!this.isSquare()) return false;
        for (let i = 1; i < this._numCols; ++i) {
            for (let j = 0; j < i; ++j) {
                if (!near(this.get(i, j), this.get(j, i), tolerance)) return false;
            }
        }
        return true;
    }
    static random(numRows: number, numCols: number): Matrix {
        let data = [];
        for (let i = 0; i < numRows * numCols; i++)
            data.push(Math.random());
        return new Matrix(data, numRows, numCols);
    }
    static empty(numRows: number, numCols: number): Matrix {
        let data: number[];
        (data = []).length = numRows * numCols;
        data.fill(0);
        return new Matrix(data, numRows, numCols);
    }
    copy(): Matrix {
        return new Matrix(this.data.slice(), this._numRows, this._numCols);
    }
    clone(): Matrix {
        return new Matrix(this.data.slice(), this._numRows, this._numCols);
    }
    width(): number {
        return this._numCols;
    }
    height(): number {
        return this._numRows;
    }
    static identity(size: number): Matrix {
        let data: number[];
        (data = []).length = size * size;
        data.fill(0);
        for (let i = 0; i < size; i++) {
            data[i + i * size] = 1;
        }
        return new Matrix(data, size, size);
    }
    static mul(a: Matrix, b: Matrix): Matrix {
        assert(a._numCols == b._numRows, "Matrix dimensions aren't compatible");
        let result = Matrix.empty(b._numCols, a._numRows);
        //for each cell in the result
        for (let j = 0; j < a._numRows; j++) {
            for (let i = 0; i < b._numCols; i++) {
                let value = 0;
                for (let k = 0; k < a._numCols; k++) {
                    value += a.get(j, k) * b.get(k, i);
                }
                result.set(j, i, value);
            }
        }
        return result;
    }
    static postMulVec(a: Matrix, b: vector): vector {
        assert(a._numCols == b.data.length, "Width of matrix isn't compatible with vector's length");
        let result = vector.empty(a._numRows);
        for (let j = 0; j < a._numRows; j++) {
            let v = 0;
            for (let i = 0; i < a._numCols; i++) {
                v += a.get(j, i) * b.get(i);
            }
            result.set(j, v);
        }
        return result;
    }
    static preMulVec(a: Matrix, b: vector): vector {
        assert(a._numRows == b.data.length, "Width of matrix isn't compatible with vector's length");
        let result = vector.empty(a._numCols);
        for (let i = 0; i < a._numCols; i++) {
            let v = 0;
            for (let j = 0; j < a._numRows; j++) {
                v += a.get(j, i) * b.get(j);
            }
            result.set(i, v);
        }
        return result;
    }
    get(row: number, column: number): number {
        assert(row < this._numRows && row >= 0, "Invalid row");
        assert(column < this._numCols && column >= 0, "Invalid column");
        return this.data[row * this._numCols + column];
    }
    set(row: number, column: number, value: number): void {
        assert(row < this._numRows && row >= 0, "Invalid row");
        assert(column < this._numCols && column >= 0, "Invalid column");
        this.data[row * this._numCols + column] = value;
    }
    transposeInPlace(): Matrix {
        assert(this.isSquare(), "Non-square matrix");
        for (let row = 0; row < this._numRows; row++) {
            for (let col = row + 1; col < this._numCols; col++) {
                swap(this.data, col + row * this._numCols, col * this._numRows + row);
            }
        }
        return this;
    }
    transpose(): Matrix {
        let result = [];
        for (let col = 0; col < this._numCols; col++) {
            for (let row = 0; row < this._numRows; row++) {
                result.push(this.data[col + row * this._numCols]);
            }
        }
        return new Matrix(result, this._numCols, this._numRows);
    }
    // partial pivot gauss elimination
    static solve(A: Matrix, b: vector): vector {
        assert(A._numCols == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A._numCols == A._numRows, "Non-square matrix");
        var rang = b.size();
        var x = vector.empty(rang);
        let epsilon = 0.001
        var indexes = new Array(rang);
        for (var i = 0; i < rang; i++) {
            indexes[i] = i;
        }
        for (var l = 0; l < rang; l++) {
            var max = l;
            for (var i = l + 1; i < rang; i++) {
                if (Math.abs(A.get(indexes[i], l)) > Math.abs(A.get(indexes[max], l)))
                    max = i;
            }
            if (max != l) {
                var temp = indexes[l];
                indexes[l] = indexes[max];
                indexes[max] = temp;
            }
            if (Math.abs(A.get(indexes[l], l)) < epsilon) {
                for (var i = 0; i < rang; i++)
                    x.set(i, 0.0);
                return x;
            }
            for (var i = l + 1; i < rang; i++)
                A.set(indexes[l], i, A.get(indexes[l], i) / A.get(indexes[l], l));
            b.set(indexes[l], b.get(indexes[l]) / A.get(indexes[l], l));
            A.set(indexes[l], l, 1);

            for (var i = l + 1; i < rang; i++) {
                for (var k = l + 1; k < rang; k++)
                    A.set(indexes[i], k, A.get(indexes[i], k) - A.get(indexes[i], l) * A.get(indexes[l], k));
                b.set(indexes[i], b.get(indexes[i]) - A.get(indexes[i], l) * b.get(indexes[l]));
                A.set(indexes[i], l, 0);
            }
        }
        x.set(rang - 1, b.get(indexes[rang - 1]));
        for (var i = rang - 2; i > -1; i--) {
            var k = 0.;
            for (var j = i + 1; j < rang; j++) {
                k = k + A.get(indexes[i], j) * x.get(j);
            }
            x.set(i, b.get(indexes[i]) - k);
        }
        return x;
    }
    static fromTriplets(triplets: Triplet[], numRows: number, numCols: number): Matrix {
        let result = Matrix.empty(numRows, numCols);
        for (const triplet of triplets)
            result.set(triplet.row, triplet.column, triplet.value);
        return result;
    }
    pseudoInverse(): Matrix {
        throw new Error("Not implemented");
    }
    inverse(tolerance: number = SmallTolerance): Matrix {
        return FullPivLU.solveMatrix(this, Matrix.identity(this.width()), tolerance);
    }
    // analytic solution
    inverseNaive(): Matrix {
        assert(this.isSquare(), "Non-square matrix");
        assert(this.width() < 7, "Factorial complexity");
        let d = this.determinantNaive();
        let result = Matrix.empty(this._numCols, this._numRows);
        let minors = this.clone();
        for (let column = 0; column < this._numCols; ++column) {
            for (let row = 0; row < this._numRows; ++row) {
                let minorValue = this.minor(column, row);
                minors.set(row, column, ((row + column) & 1) ? -minorValue : minorValue);
                minorValue = minorValue / d;
                result.set(row, column, (row + column) & 1 ? -minorValue : minorValue);
            }
        }
        console.log(`Minor ${minors.toString()} / ${d}`);
        return result;
        /*assert(this._numCols == this._numRows, "Non-square matrix");
        let result = this.copy();
        for (let i = 0; i < this._numCols; i++) {
            let v = vector.empty(this._numCols);
            v.set(i, 1);
            let column = Matrix.solve(this.copy(), v);
            for (let j = 0; j < this._numRows; j++) {
                result.set(j, i, column.get(j));
            }
        }
        return result;*/
    }
    determinantNaive() {
        assert(this.isSquare(), "Non-square matrix");
        assert(this.width() < 7, "Factorial complexity");
        if (this._numCols != this._numRows) return 0.0;

        //console.log("Determinant computation");
        //console.log(`Matrix: ${this.toString()}`);

        const self = this;
        const step = (r: number, c: number, rowList: number[], colList: number[], depth: number): number => {
            const curValue = self.get(rowList[r], colList[c]);
            if (Math.abs(curValue) < SmallestTolerance) return 0.0;

            let result = 0.0;
            if (self._numCols - depth > 0) {
                swap(rowList, r, self._numCols - depth);
                swap(colList, c, self._numCols - depth);
                for (let i = 0; i < self._numCols - depth; ++i)
                    result += step(0, i, rowList, colList, depth + 1);
                swap(rowList, r, self._numCols - depth);
                swap(colList, c, self._numCols - depth);
            } else {
                return curValue;
            }
            const curNumCols = this._numCols - depth;
            if (Math.max(curNumCols - r - 1, 0) + Math.max(curNumCols - c - 1, 0) & 1)
                result = -result;

            /*if (depth == 1) {
                console.log(`Step (${depth}) [row:  ${rowList[r]}, col: ${colList[c]}]`);
                swap(rowList, r, self._numCols - depth);
                swap(colList, c, self._numCols - depth);
                let subMatrix = Matrix.empty(self._numCols - depth, self._numRows - depth);
                for (let i = 0; i < self._numCols - depth; ++i) {
                    for (let j = 0; j < self._numRows - depth; ++j) {
                        const curValue = self.get(rowList[j], colList[i]);
                        subMatrix.set(j, i, curValue);
                    }
                }
                swap(rowList, r, self._numCols - depth);
                swap(colList, c, self._numCols - depth);
                console.log(`SubMatrix ${subMatrix.toString()}`);
                console.log(`PreResult ${result}`);
            }*/
            result *= curValue * (((r + c) & 1) == 1 ? -1 : 1);
            /*if (depth == 1)
                console.log(`Result ${result}`);*/

            return result;
        };
        let rowList = new Array(this._numCols);
        let colList = new Array(this._numCols);
        for (let i = 0; i < this._numCols; ++i)
            rowList[i] = colList[i] = i;
        let result = 0.0;
        for (let i = 0; i < this._numCols; ++i)
            result += step(0, i, rowList, colList, 1);
        return result;
    }
    minor(r: number, c: number): number {
        assert(this.isSquare(), "Non-square matrix");
        if (this._numCols == 1) return 1.0;
        let subMatrix = Matrix.empty(this._numCols - 1, this._numRows - 1);
        for (let row = 0, i = 0; row < this._numRows; ++row) {
            if (row == r) continue;
            for (let col = 0, j = 0; col < this._numCols; ++col) {
                if (col == c) continue;
                subMatrix.set(i, j, this.get(row, col));
                ++j;
            }
            ++i;
        }
        return subMatrix.determinantNaive();
    }
    print(fractionDigits: number): string {
        if (!fractionDigits)
            fractionDigits = 4;
        let result = "";
        for (let j = 0; j < this._numRows; j++) {
            result += j > 0 ? "\n| " : "| ";
            for (let i = 0; i < this._numCols; i++)
                result += this.data[i + j * this._numCols].toFixed(fractionDigits) + " "
            result += "|";
        }
        return result;
    }
    toString(): string {
        let result = "[";
        for (let j = 0; j < this._numRows; ++j) {
            if (j != 0)
                result += ",";
            result += "\n\t[";
            for (let i = 0; i < this._numCols; ++i) {
                if (i != 0) {
                    result += ", "
                }
                result += this.get(j, i).toFixed(4);
            }
            result += "]";
        }
        return result + "\n]";
    }
    static near(a: mat, b: mat, threshold: number = SmallTolerance): boolean {
        assert(a.numRows() == b.numRows() && a.numCols() == b.numCols(), "Sizes don't match");
        for (let i = 0; i < a.data.length; i++) {
            if (Math.abs(a.data[i] - b.data[i]) > threshold)
                return false;
        }
        return true;
    }
    // Frobenius norm
    static l2Distance(a: mat, b: mat): number {
        assert(a.numRows() == b.numRows() && a.numCols() == b.numCols(), "Sizes don't match");
        let result = 0.0;
        for (let i = 0; i < a.data.length; ++i) {
            const value = a.data[i] - b.data[i];
            result = value * value;
        }
        return result
    }
    // max norm
    static lInfDistance(a: mat, b: mat): number {
        assert(a.numRows() == b.numRows() && a.numCols() == b.numCols(), "Sizes don't match");
        let result = 0.0;
        for (let i = 0; i < a.data.length; ++i) {
            const value = a.data[i] - b.data[i];
            result = Math.max(result, Math.abs(value));
        }
        return result
    }
    getColumn(col: number): vector {
        let values = vector.empty(this.numRows());
        for (let i = 0; i < this.numRows(); ++i)
            values.set(i, this.get(i, col));
        return values;
    }
    getRow(row: number): vector {
        let values = vector.empty(this.numCols());
        for (let i = 0; i < this.numRows(); ++i)
            values.set(i, this.get(row, i));
        return values;
    }
    setColumn(col: number, values: vector): void {
        assert(values.size() == this.numRows(), "Invalid size");
        for (let i = 0; i < values.size(); ++i)
            this.set(i, col, values.get(i));
    }
    setRow(row: number, values: vector): void {
        assert(values.size() == this.numCols(), "Invalid size");
        for (let i = 0; i < values.size(); ++i)
            this.set(row, i, values.get(i));
    }
    static sub(a: Matrix, b: Matrix) {
        let result = a.clone();
        result.subSelf(b);
        return result;
    }
    static add(a: Matrix, b: Matrix) {
        let result = a.clone();
        result.addSelf(b);
        return result;
    }
    makeSquare(): Matrix {
        if (this.isSquare()) return this.clone();
        let size = Math.max(this.width(), this.height());
        let result = Matrix.empty(size, size);
        for (let j = 0; j < this.height(); ++j) {
            for (let i = 0; i < this.width(); ++i)
                result.set(j, i, this.get(j, i));
        }
        return result;
    }
}