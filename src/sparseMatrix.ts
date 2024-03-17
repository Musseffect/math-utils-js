import AbstractMatrix from "./abstractMatrix";
import Matrix from "./denseMatrix";
import { SparseVector, SparseVectorElement } from "./sparseVector";
import Triplet from "./triplet";
import { assert, SmallTolerance, SmallestTolerance } from "./utils";
import Vector from "./vector";

export class CellRef {
    private matrix: SparseMatrixCSR;
    private idx: number;
    private rowIdx: number;
    private colIdx: number;
    constructor(matrix: SparseMatrixCSR, idx: number, rowIdx: number, colIdx: number) {
        this.matrix = matrix;
        this.idx = idx;
        this.rowIdx = rowIdx;
        this.colIdx = colIdx;
    }
    outerIdx(): number {
        return this.idx;
    }
    row(): number {
        return this.rowIdx;
    }
    column(): number {
        return this.colIdx;
    }
    get(): number {
        return this.matrix.nonZeroElement(this.idx);
    }
    set(value: number): void {
        this.matrix.setNonZeroElement(this.idx, value);
    }
}
// list of lists
// dictionary of keys
// sorted triples

interface ValuePair {
    value1: number;
    value2: number;
    colIdx: number;
}

// sorted triplets, row major
export class SparseMatrixTriplets {
    private triplets: Triplet[];
    private width: number;
    private height: number;
    constructor(width: number, height: number) {
        throw new Error("Not implemented");
    }
    get(row: number, col: number) {
        throw new Error("Not implemented");
    }
    set(row: number, col: number, value: number) {
        // binary search
        throw new Error("Not implemented");
    }
}

class SparseMatrixTwinRowIterator {
    private m1: SparseMatrixCSR;
    private m2: SparseMatrixCSR;
    private rowIdx: number;
    private it1: number;
    private it2: number;
    constructor(m1: SparseMatrixCSR, m2: SparseMatrixCSR, rowIdx: number) {
        this.m1 = m1;
        this.m2 = m2;
        this.rowIdx = rowIdx;
        this.it1 = this.m1.outerStart(this.rowIdx);
        this.it2 = this.m2.outerStart(this.rowIdx);
    }
    isDone(): boolean {
        const isValidIt1 = this.it1 < this.m1.outerStart(this.rowIdx + 1);
        const isValidIt2 = this.it2 < this.m2.outerStart(this.rowIdx + 1);
        return !isValidIt1 && !isValidIt2;
    }
    advance(): ValuePair {
        const isValidIt1 = this.it1 < this.m1.outerStart(this.rowIdx + 1);
        const isValidIt2 = this.it2 < this.m2.outerStart(this.rowIdx + 1);
        let col1 = isValidIt1 ? this.m1.innerIndex(this.it1) : this.m1.width();
        let col2 = isValidIt2 ? this.m2.innerIndex(this.it2) : this.m2.width();
        let colIdx = 0;
        let value1 = 0.0;
        let value2 = 0.0;
        if (col1 <= col2) {
            colIdx = col1;
            value1 = this.m1.nonZeroElement(this.it1);
            ++this.it1;
        }
        if (col2 <= col1) {
            colIdx = col2;
            value2 = this.m1.nonZeroElement(this.it2);
            ++this.it2;
        }
        return { value1, value2, colIdx };
    }
}

// iterate values in a row
export class SparseMatrixRowIterator {
    private m: SparseMatrixCSR;
    private rowIdx: number;
    private it: number;
    constructor(m: SparseMatrixCSR, rowIdx: number) {
        this.m = m;
        this.rowIdx = rowIdx;
        this.it = this.m.outerStart(this.rowIdx);
    }
    isDone(): boolean {
        return this.it >= this.m.outerStart(this.rowIdx + 1);
    }
    advance(): { value: number, colIdx: number } {
        let colIdx = this.m.innerIndex(this.it);
        let value = this.m.nonZeroElement(this.it);
        ++this.it;
        return { value, colIdx };
    }
}
// CSR
// todo: CSC
// todo: track state of nonZeroElements array with state:number variable
export class SparseMatrixCSR extends AbstractMatrix {
    // column indices sizeof(NNZ)
    protected innerIndices: number[];
    // values sizeof(NNZ)
    protected nonZeroElements: number[];
    // row starts
    protected outerStarts: number[];

    constructor(numRows: number, numCols: number) {
        super(numRows, numCols);
        this.innerIndices = [];
        this.nonZeroElements = [];
        this.outerStarts = [];
    }
    clone(): SparseMatrixCSR {
        let result = new SparseMatrixCSR(this.numRows(), this.numCols());
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
    // TODO: add test
    static near(m1: SparseMatrixCSR, m2: SparseMatrixCSR, tolerance: number = SmallTolerance): boolean {
        assert(m1.width() == m2.width() && m1.height() == m2.height(), "Incompatible sizes");
        for (let row = 0; row < m1.numRows(); ++row) {

            const UseIterator = false;
            if (UseIterator) {
                let it = new SparseMatrixTwinRowIterator(m1, m2, row);
                while (!it.isDone()) {
                    let value = it.advance();
                    if (Math.abs(value.value1 - value.value2) > tolerance)
                        return false;
                }
            }
            else {
                let it1 = m1.outerStarts[row];
                let it2 = m2.outerStarts[row];
                let isValidIt1 = it1 < m1.outerStarts[row + 1];
                let isValidIt2 = it2 < m2.outerStarts[row + 1];
                // TODO: make some sort of iterator for this
                while (isValidIt1 || isValidIt2) {
                    let col1 = isValidIt1 ? m1.innerIndices[it1] : m1.width();
                    let col2 = isValidIt2 ? m2.innerIndices[it2] : m2.width();
                    let value1 = 0.0;
                    let value2 = 0.0;
                    if (col1 <= col2) {
                        value1 = m1.nonZeroElements[it1];
                        ++it1;
                    }
                    if (col2 <= col1) {
                        value2 = m1.nonZeroElements[it2];
                        ++it2;
                    }
                    if (Math.abs(value1 - value2) > tolerance)
                        return false;
                    isValidIt1 = it1 < m1.outerStarts[row + 1];
                    isValidIt2 = it2 < m2.outerStarts[row + 1];
                }
            }
        }
        return true;
    }
    static identity(size: number): SparseMatrixCSR {
        let result = new SparseMatrixCSR(size, size);
        for (let row = 0; row < size; ++row) {
            result.outerStarts.push(result.nonZeroElements.length);
            result.innerIndices.push(row);
            result.nonZeroElements.push(1);
        }
        result.outerStarts.push(result.nonZeroElements.length);
        return result;
    }
    static add(m1: SparseMatrixCSR, m2: SparseMatrixCSR): SparseMatrixCSR {
        assert(m1.width() == m2.width() && m1.height() == m2.height(), "Incompatible sizes");
        let result = new SparseMatrixCSR(m1.width(), m2.width());
        for (let row = 0; row < m1.numRows(); ++row) {
            result.outerStarts.push(result.nonZeroElements.length);
            let it = new SparseMatrixTwinRowIterator(m1, m2, row);
            while (!it.isDone()) {
                let value = it.advance();
                result.innerIndices.push(value.colIdx);
                result.nonZeroElements.push(value.value1 + value.value2);
            }
        }
        result.outerStarts.push(result.nonZeroElements.length);
        return result;
    }
    static sub(m1: SparseMatrixCSR, m2: SparseMatrixCSR): SparseMatrixCSR {
        assert(m1.width() == m2.width() && m1.height() == m2.height(), "Incompatible sizes");
        let result = new SparseMatrixCSR(m1.width(), m2.width());
        for (let row = 0; row < m1.numRows(); ++row) {
            result.outerStarts.push(result.nonZeroElements.length);
            let it = new SparseMatrixTwinRowIterator(m1, m2, row);
            while (!it.isDone()) {
                let value = it.advance();
                result.innerIndices.push(value.colIdx);
                result.nonZeroElements.push(value.value1 - value.value2);
            }
        }
        result.outerStarts.push(result.nonZeroElements.length);
        return result;
    }
    static mul(a: SparseMatrixCSR, b: SparseMatrixCSR): SparseMatrixCSR {
        // transpose second matrix and sum multiplication of row elements
        let bTransposed = b.transpose();
        let result = new SparseMatrixCSR(a.numRows(), b.numCols());
        for (let row = 0; row < a.numRows(); ++row) {
            result.outerStarts.push(result.nonZeroElements.length);
            const aNumElementsInRow = a.outerStarts[row + 1] - a.outerStarts[row];
            if (aNumElementsInRow == 0) continue;
            for (let col = 0; col < b.numCols(); ++col) {
                let aColIt = 0;
                let value = 0.0;
                for (let bRowIt = bTransposed.outerStarts[col]; bRowIt < bTransposed.outerStarts[col + 1]; ++bRowIt) {
                    const bRowIdx = bTransposed.innerIndices[bRowIt];
                    while (aColIt < aNumElementsInRow && a.innerIndices[aColIt] < bRowIdx)
                        ++aColIt;
                    if (aColIt == aNumElementsInRow) break;
                    if (a.innerIndices[aColIt] > bRowIdx) continue;
                    const bRowValue = bTransposed.nonZeroElements[bRowIt];
                    value += bRowValue * a.nonZeroElements[aColIt];
                }
                if (value != 0.0) {
                    result.nonZeroElements.push(value);
                    result.innerIndices.push(col);
                }
            }
        }
        result.outerStarts.push(result.nonZeroElements.length);
        return result;
    }
    static kroneckerProduct(a: SparseMatrixCSR, b: SparseMatrixCSR): SparseMatrixCSR {
        let result = new SparseMatrixCSR(a.numRows() * b.numRows(), a.numCols() * b.numCols());
        for (let aRow = 0, resRow = 0; aRow < a.numRows(); ++aRow) {
            for (let bRow = 0; bRow < b.numRows(); ++bRow, ++resRow) {
                result.outerStarts.push(result.nonZeroElements.length);
                for (let aOuterIdxIt = a.outerStarts[aRow]; aOuterIdxIt < a.outerStarts[aRow + 1]; ++aOuterIdxIt) {
                    let aCol = a.innerIndices[aOuterIdxIt];
                    let aVal = a.nonZeroElements[aOuterIdxIt];
                    for (let bOuterIdxIt = b.outerStarts[bRow]; bOuterIdxIt < b.outerStarts[bRow + 1]; ++bOuterIdxIt) {
                        let bCol = b.innerIndices[bOuterIdxIt];
                        let bVal = b.nonZeroElements[bOuterIdxIt];
                        let resCol = aCol * b.numCols() + bCol;
                        result.innerIndices.push(resCol);
                        result.nonZeroElements.push(aVal * bVal);
                    }
                }
            }
        }
        result.outerStarts.push(result.nonZeroElements.length);
        return result;
    }
    static postMulSparse(m: SparseMatrixCSR, v: SparseVector): SparseVector {
        assert(v.size() == m.numCols(), "Sizes don't match");
        let result: SparseVectorElement[] = [];
        for (let rowIdx = 0; rowIdx < m.numRows(); ++rowIdx) {
            let value = 0.0;
            let curVectorIdxIt = 0;
            for (let outerIdxIt = m.outerStarts[rowIdx]; outerIdxIt < m.outerStarts[rowIdx + 1]; ++outerIdxIt) {
                let colIdx = m.innerIndices[outerIdxIt];
                while (curVectorIdxIt < v.elements.length && v.elements[curVectorIdxIt].index < colIdx)
                    ++curVectorIdxIt;
                if (curVectorIdxIt == v.elements.length || v.elements[curVectorIdxIt].index > colIdx) continue;
                value += v.elements[curVectorIdxIt].value * m.nonZeroElements[outerIdxIt];
            }
            if (value != 0.0)
                result.push({ index: rowIdx, value: value });
        }
        return new SparseVector(m.numRows(), result);
    }
    static preMulSparse(v: SparseVector, m: SparseMatrixCSR): SparseVector {
        assert(v.size() == m.numRows(), "Sizes don't match");
        let transposed = m.transpose();
        return SparseMatrixCSR.postMulSparse(transposed, v);
    }
    static postMul(m: SparseMatrixCSR, v: Vector): Vector {
        assert(v.size() == m.numCols(), "Sizes don't match");
        let result: Vector = Vector.empty(m.numRows());
        for (let rowIdx = 0; rowIdx < m.numRows(); ++rowIdx) {
            let value = 0.0;
            for (let outerIdxIt = m.outerStarts[rowIdx]; outerIdxIt < m.outerStarts[rowIdx + 1]; ++outerIdxIt) {
                value += v.get(m.innerIndices[outerIdxIt]) * m.nonZeroElements[outerIdxIt];
            }
            result.set(rowIdx, value);
        }
        return result;
    }
    static preMul(v: Vector, m: SparseMatrixCSR): Vector {
        assert(v.size() == m.numRows(), "Sizes don't match");
        let transposed = m.transpose();
        return SparseMatrixCSR.postMul(transposed, v);
    }
    scale(scalar: number): SparseMatrixCSR {
        for (let i = 0; i < this.nonZeroElements.length; ++i)
            this.nonZeroElements[i] *= scalar;
        return this;
    }
    static scale(m: SparseMatrixCSR, scalar: number): SparseMatrixCSR {
        return m.clone().scale(scalar);
    }
    rowVector(row: number): SparseVector {
        throw new Error("Method not implemented.");
    }
    columnVector(column: number): SparseVector {
        throw new Error("Method not implemented.");
    }
    transpose(): SparseMatrixCSR {
        let result = new SparseMatrixCSR(this.numCols(), this.numRows());
        result.outerStarts = Array(this.numCols() + 1).fill(0);
        result.nonZeroElements = Array(this.nonZeroElements.length);
        result.innerIndices = Array(this.nonZeroElements.length);
        // calc number of entries per column/ row
        for (let i = 0; i < this.nonZeroElements.length; ++i)
            result.outerStarts[this.innerIndices[i]]++;
        // for outerStarts
        for (let sum = 0, col = 0; col < this.numCols(); ++col) {
            let curValue = result.outerStarts[col];
            result.outerStarts[col] = sum;
            sum += curValue;
        }
        // outerStarts tracks innerIdxit for transposed elements
        for (let rowIdx = 0; rowIdx < this.numRows(); ++rowIdx) {
            for (let innerIdxIt = this.outerStarts[rowIdx]; innerIdxIt < this.outerStarts[rowIdx + 1]; ++innerIdxIt) {
                let colIdx = this.innerIndices[innerIdxIt];
                let value = this.nonZeroElements[innerIdxIt];
                let resultInnerIdxIt = result.outerStarts[colIdx];
                result.outerStarts[colIdx]++;
                result.innerIndices[resultInnerIdxIt] = rowIdx;
                result.nonZeroElements[resultInnerIdxIt] = value;
            }
        }
        // restore proper format for outerStarts
        for (let col = 0, prev = 0; col <= this.numCols(); ++col) {
            let temp = result.outerStarts[col];
            result.outerStarts[col] = prev;
            prev = temp;
        }
        return result;
    }
    static fromDense(dense: Matrix, tolerance: number): SparseMatrixCSR {
        throw new Error("Method not implemented.");
    }
    static fromTriplets(triplets: Triplet[], numRows: number, numCols: number, tolerance: number = SmallestTolerance): SparseMatrixCSR {
        // sorted in ascending "row by row" order
        triplets.sort((a: Triplet, b: Triplet) => {
            let rowSign = a.row - b.row;
            if (rowSign != 0) return rowSign;
            return a.column - b.column;
        });
        let result = new SparseMatrixCSR(numRows, numCols);
        if (triplets.length == 0) return result;
        result.outerStarts.push(0);
        let currentRow = triplets[0].row;
        for (let i = 0; i < triplets.length; ++i) {
            assert(triplets[i].row <= numRows && triplets[i].row >= 0, "Invalid row index");
            assert(triplets[i].column <= numCols && triplets[i].column >= 0, "Invalid column index");
            if (Math.abs(triplets[i].value) < tolerance) continue
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
        return result;
    }
    toDense(): Matrix {
        let result = Matrix.empty(this.numRows(), this.numCols());
        for (let rowIdx = 0; rowIdx < this.numRows(); ++rowIdx) {
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
    numNonZeroes(): number {
        return this.nonZeroElements.length;
    }
    coeffRef(row: number, column: number): CellRef {
        let left = this.outerStarts[row];
        let right = this.outerStarts[row + 1];
        while (left != right) {
            let middleIt = Math.floor((left + right) / 2);
            const middleCol = this.innerIndices[middleIt];
            if (middleCol < column) {
                left = middleIt + 1;
            } else if (middleCol > column) {
                right = middleIt;
            } else {
                return new CellRef(this, middleIt, row, column);
            }
        }
        return null;
    }
    erase(cell: CellRef) {
        for (let row = cell.row() + 1; row <= this.numRows(); ++row)
            this.outerStarts[row]--;
        this.innerIndices.splice(cell.outerIdx(), 1);
        this.nonZeroElements.splice(cell.outerIdx(), 1);
    }
    set(row: number, column: number, value: number): void {
        let left = this.outerStarts[row];
        let right = this.outerStarts[row + 1];
        while (left != right) {
            let middleIt = Math.floor((left + right) / 2);
            const middleCol = this.innerIndices[middleIt];
            if (middleCol < column) {
                left = middleIt + 1;
            } else if (middleCol > column) {
                right = middleIt;
            } else {
                this.nonZeroElements[middleIt] = value;
                return;
            }
        }
        if (this.outerStarts[row] == this.outerStarts[row + 1]) {
            this.innerIndices.splice(left, 0, value);
        } else {
            throw new Error("Not implemented");
        }
        for (let rowIdx = row + 1; row <= this.numRows(); ++rowIdx)
            this.outerStarts[rowIdx]++;
    }
    // TODO: test
    get(row: number, column: number): number {
        let left = this.outerStarts[row];
        let right = this.outerStarts[row + 1];
        while (left != right) {
            let middleIt = Math.floor((left + right) / 2);
            const middleCol = this.innerIndices[middleIt];
            if (middleCol < column) {
                left = middleIt + 1;
            } else if (middleCol > column) {
                right = middleIt;
            } else {
                return this.nonZeroElements[middleIt];
            }
        }
        return 0.0;
    }
    determinant(): number {
        throw new Error("Method not implemented.");
    }
    inverse(): any {
        throw new Error("Method not implemented.");
    }
    nonZeroElement(idx: number): number {
        return this.nonZeroElements[idx];
    }
    innerIndex(idx: number): number {
        return this.innerIndices[idx];
    }
    outerStart(row: number): number {
        return this.outerStarts[row];
    }
    setNonZeroElement(idx: number, value: number) {
        return this.nonZeroElements[idx] = value;
    }
    isValid(): boolean {
        let result = true;
        result = result && this.outerStarts.length == this.numRows() + 1;
        result = result && this.innerIndices.length == this.nonZeroElements.length;
        result = result && this.outerStart(this.numRows()) == this.nonZeroElements.length;
        if (!result) return false;
        for (let row = 0; row < this.numRows(); ++row) {
            let prevCol = -1;
            for (let colIt = this.outerStarts[row]; colIt < this.outerStarts[row + 1]; ++colIt) {
                const col = this.innerIndex(colIt);
                if (col <= prevCol) return false;
                prevCol = col;

            }
        }
        return true;
    }
    isSymmetric(tolerance: number = SmallTolerance): boolean {
        let transpose = this.transpose();
        return SparseMatrixCSR.near(this, transpose, tolerance);
    }
}

// TODO: sparse matrix class that combines triplets and CSR