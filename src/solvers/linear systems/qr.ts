// todo

import Matrix from "../../denseMatrix";
import { assert } from "../../utils";
import Vector from "../../vector";

const SolverName = "QR";

export default class QR {
    A: Matrix = null;
    q: Matrix;
    r: Matrix;
    _makeCompact: boolean;
    public factorize(A: Matrix | null) {
        this.q = null;
        this.r = null;
        this.A = A;
        if (this.A == null) return;
        throw new Error("Not implemented");
    }
    constructor(A: Matrix | null, makeCompact: boolean = false) {
        this._makeCompact = makeCompact;
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