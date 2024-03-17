import Matrix from "../../denseMatrix";
import Vector from "../../vector";
import PartialPivLU from "../linear systems/partialPivLU";

export abstract class LeastSquaresResiduals {
    public abstract numResiduals(): number;
    public abstract numParams(): number;
    public abstract f(p: Vector, idx: number): number;
    public abstract dfdp(p: Vector, idx: number): Vector;
    public abstract ddfdpdp(p: Vector, idx: number): Matrix;
    public error(p: Vector): number {
        let error = 0;
        for (let i = 0; i < this.numResiduals(); i++)
            error += Math.pow(this.f(p, i), 2.0);
        return error;
    }
}

export class GaussNewton {
    numIters: number;
    errAbsTol: number;
    alpha: number;
    solve(residuals: LeastSquaresResiduals, p0: Vector): { solution: Vector, error: number } {
        let p = p0.clone();
        let numParams = residuals.numParams();
        for (let i = 0; i < this.numIters; i++) {
            let rhs = Vector.empty(numParams);
            let H = Matrix.empty(numParams, numParams);
            let error = 0.0;
            for (let resIdx = 0; resIdx < residuals.numResiduals(); resIdx++) {
                let r = residuals.f(p, resIdx);
                let drdp = residuals.dfdp(p, resIdx);
                error += r * r;
                rhs.addSelf(Vector.scale(drdp, -r * this.alpha));
                H.addSelf(Vector.outer(drdp, drdp));
            }
            if (error < this.errAbsTol)
                break;
            let dp = PartialPivLU.solve(H, rhs);
            p.addSelf(dp);
        }
        let error = residuals.error(p);
        return { solution: p, error: error };
    }
}