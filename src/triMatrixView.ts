import Matrix from "./denseMatrix";
import { assert } from "./utils";


export enum TriMatrixType {
    lower = 0,
    upper = 1
}
export class TriMatrix {
    private type: TriMatrixType;
    private data: number[];
    private size: number;
    constructor(size: number, type: TriMatrixType) {
        throw new Error("Not implemented");
    }
    private fromIndex(index: number): { row: number, column: number } {
        switch (this.type) {
            case TriMatrixType.lower:
                let row = Math.floor((Math.sqrt(1 + 8 * index) - 1) / 2);
                let column = index - row * (row + 1) / 2;
            case TriMatrixType.upper:
                throw new Error("Not implemented");
        }
    }
    private toIndex(row: number, column: number): number {
        switch (this.type) {
            case TriMatrixType.lower:
                if (column < row) return -1;
                return column + row * (row + 1) / 2;
            case TriMatrixType.upper:
                if (column > row) return -1;
            //return row + (this.size - column) * ((this.size - column) + 1) / 2;
        }
    }
    get(row: number, column: number): number {
        let index = this.toIndex(row, column);
        if (index < 0) return 0.0;
        return this.data[index];
    }
    set(row: number, column: number, value: number): void {
        let index = this.toIndex(row, column);
        assert(index < 0, "Invalid index");
        this.data[index] = value;
    }
}

export enum DiagonalType {
    Unit,
    Zero,
    Existing        
};

export class TriMatrixView {
    private m: Matrix;
    private type: TriMatrixType;
    private diagType: DiagonalType;
    constructor(m: Matrix, type: TriMatrixType, diagType: DiagonalType) {
        this.m = m;
        this.type = type;
        this.diagType = diagType;
    }
    get(row: number, column: number): number {
        if (row == column)
        {
            switch(this.diagType) {
                case DiagonalType.Zero: return 0.0;
                case DiagonalType.Unit: return 1.0;
            }
        }
        switch (this.type) {
            case TriMatrixType.lower:
                if (column > row)
                    return 0.0;
            case TriMatrixType.upper:
                if (column < row)
                    return 0.0;
        }
        return this.m.get(row, column);
    }
    toMatrix(): Matrix {
        let result = Matrix.empty(this.m.numRows(), this.m.numCols());
        switch (this.type) {
            case TriMatrixType.lower:
                for (let j = 0; j < result.height(); ++j) {
                    for (let i = 0; i <= j; ++i)
                        result.set(j, i, this.get(j, i));
                }
                return result;
            case TriMatrixType.upper:
                for (let j = 0; j < result.height(); ++j) {
                    for (let i = j; i <= result.width(); ++i)
                        result.set(j, i, this.get(j, i));
                }
                return result;
        }
    }
    toTriangularMatrix(): TriMatrix {
        throw new Error("Not implemented");
    }
}