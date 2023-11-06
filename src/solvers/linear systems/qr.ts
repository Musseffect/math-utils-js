// todo

import Matrix from "../../denseMatrix";
import { assert, assertFail } from "../../utils";
import Vector from "../../vector";
import { applyGivensFromLeft, applyGivensFromRight, applyHouseholderFromLeft, applyHouseholderFromRight, calcHouseholderVectorCol, givens } from "./eigenvalues";

const SolverName = "QR";

export enum ZeroingMethod {
    Givens = 0,
    Housholder = 1
};

export class QR {
    A: Matrix = null;
    q: Matrix;
    r: Matrix;
    _zeroingMethod: ZeroingMethod;
    _makeCompact: boolean;
    private factorizeGivens() {
        // make zeroes for each column from the bottom to rowIdx - 1
        for (let j = 0; j != this.A.numCols(); ++j) {
            // | a
            // | x
            // | b  
            // | 0  /|\
            // | 0 / | \
            for (let i = this.A.numRows() - 1; i > j; --i) {
                let a = this.r.get(j, j);
                let b = this.r.get(i, j);
                let givensCoeffs = givens(a, b);
                // apply to r
                applyGivensFromLeft(this.r, givensCoeffs, i, j);
                this.r.set(j, j, givensCoeffs.r);
                this.r.set(i, j, 0.0);
                // apply to q
                applyGivensFromRight(this.q, givensCoeffs, i, j);
            }
        }
    }
    private factorizeHousholder() {
        for (let col = 0; col < this.A.numCols(); ++col) {
            if (this.A.numRows() == col + 1) continue;
            // calc housholder
            let v = calcHouseholderVectorCol(this.r, col, col);
            applyHouseholderFromLeft(v, this.r, col);

            this.r.set(col, col, v.get(0));
            for (let row = col + 1; row < this.A.numRows(); ++row)
                this.r.set(row, col, 0.0);
            // R = (Q3Q2Q1)^A = QT*A
            // QT = Q3(Q2(Q1*I))
            // Q = ((I*Q1)Q2)Q3)
            // QNT = QN because householder matrix is symmetric
            applyHouseholderFromRight(v, this.q, col);
        }
    }
    public factorize(A: Matrix | null) {
        this.q = null;
        this.r = null;
        this.A = A;
        if (this.A == null) return;

        this.r = this.A.clone();
        if (this.A.numCols() > this.A.numRows())
            this.r.transposeInPlace();
        this.q = Matrix.identity(this.A.numRows());
        switch (this._zeroingMethod) {
            case ZeroingMethod.Givens:
                this.factorizeGivens();
            case ZeroingMethod.Housholder:
                this.factorizeHousholder();
            default:
                assertFail("Invalid value");
        }
        if (this.makeCompact) {
            this.r.shrinkRows(this.A.numCols());
            this.q.shrinkCols(this.A.numCols());
        }
    }
    constructor(A: Matrix | null, zeroingMethod: ZeroingMethod = ZeroingMethod.Housholder, makeCompact: boolean = false) {
        this._makeCompact = makeCompact;
        this._zeroingMethod = zeroingMethod;
        this.factorize(A);
    }
    set makeCompact(value: boolean) {
        this._makeCompact = value;
    }
    get Q(): Matrix {
        return this.q;
    }
    get R(): Matrix {
        return this.r;
    }
    solve(x: Matrix | Vector): Matrix | Vector {
        if (this.A.numRows() > this.A.numCols()) {
            // overdetermined system
            // x = R_1^-1(Q_1^T * b), solve by back substitution
        }
        else if (this.A.numCols() > this.A.numRows()) {
            // underdetermined system
            // AT=QR
            return null;
        } else {

        }
        throw new Error("Not implemented");
    }
}

// todo: column pivoting, full pivoting