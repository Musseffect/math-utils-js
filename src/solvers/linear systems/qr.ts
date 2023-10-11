// todo

import Matrix from "../../denseMatrix";
import { assert } from "../../utils";
import Vector from "../../vector";
import { ConvergenseFailureException } from "./exceptions";

const SolverName = "QR";

class QR {
    A: Matrix;
    Q: Matrix;
    R: Matrix;
    protected decompose() {
        assert(this.A.isSquare(), "Expected square matrix");
        throw new Error("Not implemented");
    }
    constructor(A: Matrix, compact: boolean) {
        this.decompose();
    }
}