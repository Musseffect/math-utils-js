import Matrix from "../../denseMatrix";
import { assert } from "../../utils";


export default class PolarDecomposition {
    // orthogonal
    q: Matrix;
    // positive semi-definite
    s: Matrix;
    A: Matrix;
    public factorize(A: Matrix | null) {
        this.A = A;
        this.q = null;
        this.s = null;
        if (!this.A) return;
        // dynamically weighted halley
        this.s = this.A.clone();
        //https://www.cs.ucdavis.edu/~bai/Winter09/nakatsukasabaigygi09.pdf
        throw new Error("Not implemented");
    }
    constructor(A: Matrix | null = null) {
        this.factorize(A);
    }
    get Q() {
        return this.q;
    }
    get S() {
        return this.s;
    }
}