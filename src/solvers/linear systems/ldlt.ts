import Matrix from "../../denseMatrix";
import { TriMatrixView } from "../../triMatrixView";
import Vector from "../../vector";


export default class LDLT {
    LDLT: Matrix;
    constructor(tolerance: number) {
    }
    factorize(A:Matrix) {
        throw new Error("Not implemented");
    }
    solve(rhs:Vector):Vector {
        throw new Error("Not implemented");
    }
    L(): TriMatrixView {
        throw new Error("Not implemented");
    }
    D() {
        throw new Error("Not implemented");
    }
}