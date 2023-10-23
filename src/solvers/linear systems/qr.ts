// todo

import Matrix from "../../denseMatrix";
import { assert } from "../../utils";

const SolverName = "QR";

export default class QR {
    A: Matrix = null;
    q: Matrix;
    r: Matrix;
    _makeCompact: boolean;
    public factorize(A: Matrix | null) {
        this.q = null;
        this.r = null;
        if (A == null) return;
        assert(this.A.isSquare(), "Expected square matrix");
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
}