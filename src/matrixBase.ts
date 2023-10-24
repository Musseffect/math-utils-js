enum TestOrder {
    RowMajor,
    ColumnMajor
}

abstract class MatrixOrder {
    abstract eval(): boolean;
    static eval2() {
        return true;
    }
};
export abstract class RowMajor extends MatrixOrder {
    static eval() {
        return true;
    }
};
export abstract class ColumnMajor extends MatrixOrder {
    static eval() {
        return false;
    }
};

class TrueType {
    static eval(): boolean {
        return true;
    }
}
class FalseType {
    static eval(): boolean {
        return false;
    }
}

type IsRowMajor<T extends RowMajor | ColumnMajor> = (T extends RowMajor ? TrueType : FalseType);
type IsRowMajor2<T extends RowMajor | ColumnMajor> = (T extends RowMajor ? true : false);

function value<T extends RowMajor | ColumnMajor>(): IsRowMajor2<T> {
    return;
}

/*abstract*/ class MatrixStorage<T extends RowMajor | ColumnMajor> {
    isRowMajor: boolean;
    constructor() {
        this.isRowMajor = value<T>();
    }
    /*abstract get(row: number, col: number): number;
    abstract set(row: number, col: number, value: number): void;
    abstract numRows(): number;
    abstract numCols(): number;*/
    toString(): boolean {
        let p: T;
        if (p instanceof RowMajor)
            return true;
        else
            return false;
    }
};

class TriangularMatrixStorage {

}

class SymmetricMatrixStorage {
}

// todo: row major, col major, triangular, diagonal, tridiagonal

abstract class MatrixBase<T extends MatrixStorage> {
    private storage: T;
    get(row: number, col: number): number {
        return this.storage.get(row, col);
    }
    set(row: number, col: number, value: number): void {
        this.storage.set(row, col, value);
    }
    numRows(): number {
        return this.storage.numRows();
    }
    numCols(): number {
        return this.storage.numCols()
    }
    width(): number {
        return this.storage.numCols();
    }
    height(): number {
        return this.storage.numRows();
    }
    isSquare(): boolean {
        return this.width() == this.height();
    }
};