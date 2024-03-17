import Matrix from "../../denseMatrix";
import { assert } from "../../utils";
import Vector from "../../vector";
import { ConvergenseFailureException } from "./exceptions";
import { makeHessenbergInplace, makeTridiagonalInplace } from "./hessenbergMatrix";

// todo: Givens rotations, complex eigenvalues
/**
 * Compute real parts of eigenvalues with Francis double step QR algorithm for arbitary square matrix
 * @param A Input matrix
 * @param numIters number of QR iterations
 * @param tolerance tolerance for zeroed elements
 * @returns array of real eigenvalues
 */
export function calcEigenvalues(A: Matrix, numIters: number, tolerance: number): number[] {
    assert(A.isSquare(), "Expected square matrix");
    let eigenvalues: number[] = new Array(A.numCols());
    if (A.numCols() == 1) {
        eigenvalues[0] = A.get(0, 0);
        return eigenvalues;
    }
    else if (A.numCols() == 2) {
        let b = A.get(0, 0) + A.get(1, 1);
        let c = A.get(0, 0) * A.get(1, 1) - A.get(0, 1) * A.get(1, 0);
        let D = (b * b - 4.0 * c);
        eigenvalues[0] = b * 0.5;
        eigenvalues[1] = b * 0.5;
        if (D > 0) {
            D = Math.sqrt(D);
            eigenvalues[0] += D * 0.5;
            eigenvalues[1] -= D * 0.5;
        }
        return eigenvalues;
    }
    A = A.clone();

    let u = Vector.empty(A.numCols());
    if (!A.isHessenberg()) {
        // turn into hessenberg matrix
        for (let i = 0; i < u.size() - 2; i++) {
            let xNorm = 0.0;
            for (let k = 0, j = i + 1; j < u.size(); j++, k++) {
                u.set(k, A.get(j, i));
                xNorm += u.get(k) * u.get(k);
            }
            let ro = -Math.sign(A.get(i + 1, i));
            let uNorm = xNorm - A.get(i + 1, i) * A.get(i + 1, i);
            u.set(0, u.get(0) - ro * Math.sqrt(xNorm));
            uNorm += u.get(0) * u.get(0);
            uNorm = Math.sqrt(uNorm);
            u.scaleSelf(1.0 / uNorm);
            let u_a = Vector.empty(u.size() - i); //uk* Ak+1:n,k:n
            // premultiply by Q_i
            for (let j = i; j < u.size(); j++) {
                let value = 0.0;
                for (let k = i + 1; k < u.size(); k++)
                    value += u.get(k - i - 1) * A.get(k, j);
                u_a.set(j - i, value);
            }

            for (let j = i + 1; j < u.size(); j++) {
                for (let k = i; k < u.size(); k++)
                    A.set(j, k, A.get(j, k) - u.get(j - i - 1) * 2.0 * u_a.get(k - i));
            }
            // postmultiply by Q_i
            u_a = Vector.empty(u.size());
            for (let j = 0; j < u.size(); j++) {
                let value = 0.0;
                for (let k = i + 1; k < u.size(); k++)
                    value += u.get(k - i - 1) * A.get(j, k);
                u_a.set(j, value);
            }

            for (let j = 0; j < u.size(); j++) {
                for (let k = i + 1; k < u.size(); k++)
                    A.set(j, k, A.get(j, k) - 2.0 * u_a.get(j) * u.get(k - i - 1));
            }
        }
    }
    //console.log(`Hessenberg ${A.toString()}`);
    // A = makeHessenberg(A);

    for (let i = 0; i < u.size() - 2; i++) {
        for (let j = i + 2; j < u.size(); j++)
            A.set(j, i, 0.0);
    }
    // Francis double step QR
    let iter = 0;
    for (let p = u.size() - 1; p > 1;) {
        let q = p - 1;
        let s = A.get(q, q) + A.get(p, p);
        let t = A.get(q, q) * A.get(p, p) - A.get(p, q) * A.get(q, p);
        let x = A.get(0, 0) * A.get(0, 0) + A.get(0, 1) * A.get(1, 0) - s * A.get(0, 0) + t;
        let y = A.get(1, 0) * (A.get(0, 0) + A.get(1, 1) - s);
        let z = A.get(1, 0) * A.get(2, 1);
        for (let k = 0; k <= p - 2; k++) {
            let r = Math.max(0, k - 1);
            let p_v = new Vector([x, y, z]);
            let ro = -Math.sign(x);
            p_v.set(0, p_v.get(0) - ro * p_v.l2Norm());
            p_v.normalize();

            let p_t = new Array(u.size() - r);
            for (let j = r, m = 0; j < u.size(); j++, m++) {
                let temp = 0.0;
                for (let i = k, l = 0; l < 3; i++, l++)
                    temp += p_v.get(l) * A.get(i, j);
                p_t[m] = temp;
            }
            for (let j = k, l = 0; l < 3; j++, l++) {
                for (let i = r; i < u.size(); i++)
                    A.set(j, i, A.get(j, i) - 2.0 * p_v.get(l) * p_t[i - r]);
            }
            r = Math.min(k + 3, p);
            p_t = new Array(r + 1);
            for (let j = 0; j <= r; j++) {
                let value = 0.0;
                for (let i = k, l = 0; l < 3; i++, l++)
                    value += p_v.get(l) * A.get(j, i);
                p_t[j] = value;
            }

            for (let i = 0; i <= r; i++) {
                for (let j = k, l = 0; l < 3; j++, l++)
                    A.set(i, j, A.get(i, j) - 2.0 * p_v.get(l) * p_t[i]);
            }
            x = A.get(k + 1, k);
            y = A.get(k + 2, k);
            if (k < p - 2)
                z = A.get(k + 3, k);
        }

        let p_v = new Vector([x, y]);
        let ro = -Math.sign(x);
        p_v.set(0, p_v.get(0) - ro * p_v.l2Norm());
        p_v.normalize();

        let p_t = new Array(u.size() - p + 2);
        for (let j = p - 2, m = 0; j < u.size(); j++, m++) {
            let temp = 0.0;
            for (let i = q; i <= p; i++)
                temp += p_v.get(i - q) * A.get(i, j);
            p_t[m] = temp;
        }
        for (let i = q; i <= p; i++) {
            for (let j = p - 2, m = 0; j < u.size(); j++, m++)
                A.set(i, j, A.get(i, j) - 2.0 * p_v.get(i - q) * p_t[m]);
        }


        p_t = new Array(p + 1);
        for (let j = 0; j <= p; j++) {
            let value = 0.0;
            for (let i = p - 1, m = 0; i <= p; i++, m++)
                value += p_v.get(m) * A.get(j, i);
            p_t[j] = value;
        }

        for (let i = 0; i <= p; i++) {
            for (let j = p - 1, m = 0; j <= p; j++, m++)
                A.set(i, j, A.get(i, j) - 2.0 * p_v.get(m) * p_t[i]);
        }
        if (Math.abs(A.get(p, q)) < tolerance * (Math.abs(A.get(q, q)) + Math.abs(A.get(p, p)))) {
            A.set(p, q, 0);
            p = p - 1;
            q = p - 1;
        } else if (Math.abs(A.get(p - 1, q - 1)) < tolerance * (Math.abs(A.get(q - 1, q - 1)) + Math.abs(A.get(q, q)))) {
            A.set(p - 1, q - 1, 0);
            p = p - 2;
            q = p - 1;
        }
        iter++;
        if (iter > numIters)
            throw new ConvergenseFailureException("EignevaluesSolver");
    }
    for (let i = 0; i < A.numCols(); i++) {
        if (i > 0 && Math.abs(A.get(i, i - 1)) > tolerance * 10.0) //complex eigenvalues
        {
            let b = A.get(i - 1, i - 1) + A.get(i, i);
            let c = A.get(i - 1, i - 1) * A.get(i, i) - A.get(i - 1, i) * A.get(i, i - 1);
            let D = b * b - 4.0 * c;
            let x1 = b * 0.5;
            let x2 = b * 0.5;
            if (D > 0) {
                D = Math.sqrt(D);
                x1 += D * 0.5;
                x2 -= D * 0.5;
            }
            eigenvalues[i - 1] = x1;
            eigenvalues[i] = x2;
        } else {
            eigenvalues[i] = A.get(i, i);
        }
    }
    return eigenvalues;
}

// todo: use for regularization of Newton optimization method by clamping negative eigenvalues
/**
 * Q - orthogonal matrix, D - diagonal matrix of eigenvalues of A
 */

class SymmetricEigendecomposition {
    private q: Matrix = null;
    private d: Vector = null;
    private A: Matrix = null;
    constructor(A: Matrix | null = null) {
        if (A == null) return;
        this.factorize(A);
    }
    public factorize(A: Matrix) {
        assert(A.isSymmetric(), "Expected symmetric matrix");
        this.A = A;
        this.q = null;
        this.d = null;
        let T = A.clone();
        this.q = Matrix.identity(A.numCols());
        if (!T.isTridiagonal())
            makeTridiagonalInplace(T, this.q);
        // symmetric matrices have only real eigenvalues so single shift algorithm can be used
        throw new Error("Not implemented");
    }
    public get D(): Vector {
        return this.d;
    }
    public get Q(): Matrix {
        return this.q;
    }
}

class Eigendecomposition {
    private q: Matrix = null;
    private d: Matrix = null;
    private A: Matrix = null;
    constructor(A: Matrix | null = null) {
        if (A == null) return;
        this.factorize(A);
    }
    public factorize(A: Matrix) {
        this.A = A;
        this.q = null;
        this.d = A.clone();
        this.q = Matrix.identity(A.numCols());
        if (!this.d.isHessenberg(true))
            makeHessenbergInplace(this.d, this.q);
        // symmetric matrices have only real eigenvalues so single shift algorithm can be used
        throw new Error("Not implemented");
    }
    public get D(): Matrix {
        return this.d;
    }
    public get Q(): Matrix {
        return this.q;
    }
}
export function calcEigendecomposition(A: Matrix): { Q: Matrix, D: Matrix } {
    throw new Error("Not implemented");
}

// todo: Jacobi method for eigenvalues of symmetric matrix, also for singular values too