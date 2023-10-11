import Matrix from "../../denseMatrix";
import { TriMatrixView } from "../../triMatrixView";


class LDL {
    LDL: Matrix;
    A: Matrix;
    private decompose() {
        throw new Error("Not implemented'");
    }
    constructor(A: Matrix, tolerance: number) {
        this.A = A;
        this.decompose();
    }
    L(): TriMatrixView {
        throw new Error("Not implemented");
    }
    D() {
        throw new Error("Not implemented");
    }
}