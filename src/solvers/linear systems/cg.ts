import Matrix from "../../denseMatrix";
import SparseMatrix from "../../sparseMatrix";
import { assert, SmallTolerance } from "../../utils";
import Vector from "../../vector";
import { ConvergenseFailureException } from "./exceptions";

const SolverName = "'Cholesky'";

// for sparce systems
export default class CG {
    static solve(A: Matrix, b: Vector, maxIterations: number = 20, tolerance: number = SmallTolerance): Vector {
        assert(A.width() == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.isSquare(), "Non-square matrix");
        let x = Vector.empty(b.size());
        let r = Vector.sub(b, Matrix.postMulVec(A, x));
        let rDot = Vector.dot(r, r);
        if (r.lInfNorm() < tolerance) return x;
        let p = r.clone();
        let iter = 0;
        while (iter < maxIterations) {
            let ap = Matrix.postMulVec(A, p);
            let alpha = rDot / Vector.dot(p, ap);
            x.addSelf(Vector.scale(p, alpha));
            r.subSelf(Vector.scale(ap, alpha));
            let rDotNew = Vector.dot(r, r);
            if (r.lInfNorm() < tolerance) return x;
            let beta = rDotNew / rDot;
            p = p.scaleSelf(beta).addSelf(r);
            rDot = rDotNew;
            ++iter;
        }
        throw new ConvergenseFailureException(SolverName);
    }
    /*
    static solveSparse(A: SparseMatrix, b:Vector, maxIterations:number = 20, tolerance: number = SmallTolerance):Vector
    {
        assert(A.width() == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.isSquare(), "Non-square matrix");

        let x = Vector.empty(b.size());
        let r = Vector.sub(b, SparseMatrix.postMulVec(A, x));
        let rDot = Vector.dot(r, r);
        if (r.lInfNorm() < tolerance) return x;
        let p = r.clone();
        let iter = 0;
        while (iter < maxIterations) {
            let ap = SparseMatrix.postMulVec(A, p);
            let alpha = rDot / Vector.dot(p, ap);
            x.addSelf(Vector.scale(p, alpha));
            r.subSelf(Vector.scale(ap, alpha));
            let rDotNew = Vector.dot(r, r);
            if (r.lInfNorm() < tolerance) return x;
            let beta = rDotNew / rDot;
            p = p.scaleSelf(beta).addSelf(r);
            rDot = rDotNew;
            ++iter;
        }
        throw new ConvergenseFailureException(SolverName);
    }
    */
}