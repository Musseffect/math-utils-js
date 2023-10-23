

abstract class MatrixStorage {
    abstract get(row:number, col:number):number;
    abstract set(row:number, col:number, value:number):void;
    abstract numRows():number;
    abstract numCols():number;
};

// todo: row major, col major, triangular, diagonal, tridiagonal

abstract class MatrixBase<T extends MatrixStorage> {
    private storage: T;
    get(row:number, col:number): number {
        return this.storage.get(row, col);
    }
    set(row:number, col:number, value:number):void {
        this.storage.set(row, col, value);
    }
    numRows():number {
        return this.storage.numRows();
    }
    numCols():number {
        return this.storage.numCols()
    }
    width():number {
        return this.storage.numCols();
    }
    height():number {
        return this.storage.numRows();
    }
    isSquare(): boolean {
        return this.width() == this.height();
    }
};