import Matrix from "../../dense/denseMatrix";
import { SmallestTolerance, assert, assertFail, sign } from "../../utils";
import Vector from "../../dense/vector";
import { givens, applyGivensFromLeft, applyGivensFromRight, applyTransposeGivensFromRight, applyTransposeGivensFromLeft } from "./givensRotation";
import { applyHouseholderFromLeft, applyHouseholderFromRight, calcHouseholderVectorCol, calcHouseholderVectorRow } from "./hausholderReflection";
import { PermutationMatrix, PermutationType } from "../../permutationMatrix";

export enum ZeroingMethod {
    Givens = 0,
    Housholder = 1
};

export enum OrthogonalDecompositionType {
    QR,
    RQ,
    LQ,
    QL
}

export class OrthogonalDecompositionParams {
    zeroingMethod: ZeroingMethod = ZeroingMethod.Housholder;
    makeCompact: boolean = true;
    enablePivoting: boolean = true;
    normSqrTolerance: number = SmallestTolerance;
    setMethod(method: ZeroingMethod): OrthogonalDecompositionParams {
        this.zeroingMethod = method;
        return this;
    }
    setCompact(value: boolean): OrthogonalDecompositionParams {
        this.makeCompact = value;
        return this;
    }
    setPivoting(value: boolean): OrthogonalDecompositionParams {
        this.enablePivoting = value;
        return this;
    }
    setNormSqrTolerance(value: number): OrthogonalDecompositionParams {
        this.normSqrTolerance = value;
        return this;
    }
}

export class OrthogonalDecomposition {
    protected q: Matrix = null;
    protected t: Matrix = null;
    protected p: PermutationMatrix = null;
    protected params: OrthogonalDecompositionParams;
    protected determinantSign: number;
    protected type: OrthogonalDecompositionType;
    protected rank: number;

    constructor(A: Matrix | null, type: OrthogonalDecompositionType = OrthogonalDecompositionType.QR, params: OrthogonalDecompositionParams = new OrthogonalDecompositionParams()) {
        this.params = params;
        this.type = type;
        this.factorize(A);
    }
    setType(type: OrthogonalDecompositionType) {
        this.type = type;
    }
    get Q(): Matrix {
        return this.q;
    }
    get T(): Matrix {
        return this.t;
    }
    get P(): PermutationMatrix {
        return this.p;
    }
    get Rank(): number {
        return this.rank;
    }
    private makeQR(A: Matrix) {
        let r: Matrix = A.clone();
        if (A.numCols() > A.numRows())
            r.transposeInPlace();
        let q: Matrix = Matrix.identity(r.numRows());
        let p: PermutationMatrix = PermutationMatrix.identity(r.numCols(), PermutationType.Col);
        if (this.params.zeroingMethod == ZeroingMethod.Housholder) {
            let numOps = 0;
            for (let step = 0; step < r.numCols(); ++step) {
                // numRows cannot be less than numCols, when numRows == numCols last column should be skipped
                if (A.numRows() == step + 1) break;
                numOps++;
                let pivotIdx = 0;
                let normSqr = 0.0;
                for (let colIdx = step; colIdx < r.numCols(); ++colIdx) {
                    let col = p.at(colIdx);
                    let curNormSqr = 0;
                    for (let row = step; row < r.numRows(); ++row)
                        curNormSqr += Math.pow(r.get(row, col), 2);
                    if (curNormSqr > normSqr) {
                        pivotIdx = colIdx;
                        normSqr = curNormSqr;
                    }
                }
                if (normSqr <= this.params.normSqrTolerance) {
                    return;
                    // todo: check what to do in this situation
                    throw new Error("Not implemented");
                }
                r.swapColumns(step, pivotIdx);
                p.swap(step, pivotIdx);

                let v = r.subColumn(step, step, r.numRows() - step);
                let ro = -sign(v.get(0));
                if (normSqr == 0.0) return v;
                let norm = Math.sqrt(normSqr);
                let firstElement = v.get(0);
                v.set(0, v.get(0) - ro * norm);
                v.scaleSelf(1.0 / Math.sqrt(normSqr - firstElement * firstElement + v.get(0) * v.get(0)));

                applyHouseholderFromLeft(v, r, step);
                applyHouseholderFromRight(v, q, step);
            }
            this.determinantSign = (numOps & 1 ? -1.0 : 1.0);

        } else if (this.params.zeroingMethod == ZeroingMethod.Givens) {
            for (let j = 0; j != r.numCols(); ++j) {
                for (let i = r.numRows() - 1; i > j; --i) {
                    let a = r.get(j, j);
                    let b = r.get(i, j);
                    let givensCoeffs = givens(a, b);
                    applyGivensFromLeft(r, givensCoeffs, i, j);
                    applyTransposeGivensFromRight(q, givensCoeffs, i, j);
                }
            }
            this.determinantSign = 1.0;
        }
        if (this.params.makeCompact) {
            r.shrinkRows(r.numCols());
            q.shrinkCols(r.numCols());
        }
        this.q = q;
        this.t = r;
        this.p = p;
    }
    // todo: check on rectangular matrices
    // todo: check on square and rectangular matrices of different ranks
    private makeRQ(A: Matrix) {
        let r: Matrix = A.clone();
        if (A.numRows() > A.numCols())
            r.transposeInPlace();
        let q: Matrix = Matrix.identity(r.numCols());
        if (this.params.zeroingMethod == ZeroingMethod.Housholder) {
            let numOps = 0;
            for (let row = A.numRows() - 1; row > 0; --row) {
                let xNorm = 0.0;
                for (let col = 0; col <= row; ++col)
                    xNorm += Math.pow(r.get(row, col), 2);
                numOps++;
                xNorm = -sign(r.get(row, row)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorRow(r, row, 0, row + 1, row);
                applyHouseholderFromRight(v, r, 0);
                applyHouseholderFromLeft(v, q, 0);
            }
            this.determinantSign = (numOps & 1 ? -1.0 : 1.0);
        } else if (this.params.zeroingMethod == ZeroingMethod.Givens) {
            for (let row = r.numRows() - 1; row > 0; --row) {
                for (let col = 0; col < row; ++col) {
                    let i = row;
                    let j = col;
                    let givensCoeffs = givens(r.get(row, i), r.get(row, j));
                    applyGivensFromRight(r, givensCoeffs, i, j);
                    applyTransposeGivensFromLeft(q, givensCoeffs, i, j);
                }
            }
            this.determinantSign = 1.0;
        }
        if (this.params.makeCompact) {
            r.shrinkCols(r.numRows());
            q.shrinkRows(r.numRows());
        }
        this.q = q;
        this.t = r;
    }
    private makeQL(A: Matrix) {
        let q: Matrix;
        let l: Matrix;
        l = A.clone();
        q = Matrix.identity(A.numRows());
        if (this.params.zeroingMethod == ZeroingMethod.Housholder) {
            let numOps = 0;
            for (let col = l.numCols() - 1; col > 0; --col) {
                let xNorm = 0.0;
                for (let row = 0; row <= col; ++row)
                    xNorm += Math.pow(l.get(row, col), 2);
                xNorm = -sign(l.get(col, col)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorCol(l, 0, col, col + 1, col);
                applyHouseholderFromLeft(v, l, 0);
                applyHouseholderFromRight(v, q, 0);
            }
            this.determinantSign = (numOps & 1 ? -1.0 : 1.0);
        } else if (this.params.zeroingMethod == ZeroingMethod.Givens) {
            for (let col = l.numCols() - 1; col > 0; --col) {
                for (let row = 0; row < col; ++row) {
                    let i = col;
                    let j = row;
                    let givensCoeffs = givens(l.get(i, col), l.get(j, col));
                    applyTransposeGivensFromLeft(l, givensCoeffs, i, j);
                    applyGivensFromRight(q, givensCoeffs, i, j);
                }
            }
            this.determinantSign = 1.0;
        }
        this.q = q;
        this.t = l;
    }
    private makeLQ(A: Matrix) {
        let q: Matrix;
        let l: Matrix;
        l = A.clone();
        q = Matrix.identity(A.numRows());
        if (this.params.zeroingMethod == ZeroingMethod.Housholder) {
            let numOps = 0;
            for (let row = 0; row < l.numRows(); ++row) {
                let xNorm = 0.0;
                for (let col = row; col < l.numCols(); ++col)
                    xNorm += Math.pow(l.get(row, col), 2);
                xNorm = -sign(l.get(row, row)) * Math.sqrt(xNorm);
                let v = calcHouseholderVectorRow(l, row, row);
                applyHouseholderFromRight(v, l, row);
                applyHouseholderFromLeft(v, q, row);
                expect(l.get(row, row)).toBeCloseTo(xNorm);
                for (let col = row + 1; col < l.numCols(); ++col)
                    expect(l.get(row, col)).toBeCloseTo(0);
            }
            this.determinantSign = (numOps & 1 ? -1.0 : 1.0);
        } else if (this.params.zeroingMethod == ZeroingMethod.Givens) {
            for (let row = 0; row < l.numRows(); ++row) {
                for (let col = l.numCols() - 1; col > row; --col) {
                    let i = col;
                    let j = row;
                    let givensCoeffs = givens(l.get(row, j), l.get(row, i));
                    applyTransposeGivensFromRight(l, givensCoeffs, i, j);
                    applyGivensFromLeft(q, givensCoeffs, i, j);
                }
            }
            this.determinantSign = 1.0;
        }
        this.q = q;
        this.t = l;
    }
    factorize(A: Matrix) {
        if (A == null) return;
        switch (this.type) {
            case OrthogonalDecompositionType.LQ:
                /*let decomposition = new LQ(A, this.params);
                this.q = decomposition.Q;
                this.t = decomposition.L;
                this.p = decomposition.P;
                this.rank = decomposition.Rank;
                this.determinantSign = decomposition.determinantSign();*/
                this.makeLQ(A);
                break;
            case OrthogonalDecompositionType.QL:
                /*let decomposition = new QL(A, this.params);
                this.q = decomposition.Q;
                this.t = decomposition.L;
                this.p = decomposition.P;
                this.rank = decomposition.Rank;
                this.determinantSign = decomposition.determinantSign();*/
                this.makeQL(A);
                break;
            case OrthogonalDecompositionType.QR:
                let decomposition = new QR(A, this.params);
                this.q = decomposition.Q;
                this.t = decomposition.R;
                this.p = decomposition.P;
                this.rank = decomposition.Rank;
                this.determinantSign = decomposition.determinantSign();
                // this.makeQR(A); 
                break;
            case OrthogonalDecompositionType.RQ:
                /*let decomposition = new RQ(A, this.params);
                this.q = decomposition.Q;
                this.t = decomposition.R;
                this.p = decomposition.P;
                this.rank = decomposition.Rank;
                this.determinantSign = decomposition.determinantSign();*/
                this.makeRQ(A);
                break;
        }
        if (this.t != null) {
            console.log(this.t.toString());
            console.log(this.q.toString());
        }
    }
    get Params(): OrthogonalDecompositionParams {
        return this.params;
    }
    determinant(): number {
        throw new Error("Not implemented");
    }
    private solveMatrix(rhs: Matrix): Matrix {
        throw new Error("Not implemented");
    }
    private solveVector(rhs: Vector): Vector {
        throw new Error("Not implemented");
    }
    solve(rhs: Matrix | Vector): Matrix | Vector {
        if (rhs instanceof Matrix)
            return this.solveMatrix(rhs);
        return this.solveVector(rhs);
    }
    // returns inverse or pseudoinverse depending on the rank of the system
    inverse(): Matrix {
        throw new Error("Not implemented");
    }
}

export class QR {
    private A: Matrix = null;
    private q: Matrix = null;
    private r: Matrix = null;
    private rank: number = 0;
    protected p: PermutationMatrix = null;
    private qDet: number = 0;
    private params: OrthogonalDecompositionParams;

    get Params(): OrthogonalDecompositionParams {
        return this.params;
    }
    get Rank(): number {
        return this.rank;
    }
    private factorizeGivens() {
        let r: Matrix = this.r;
        this.r = null;
        let q: Matrix = Matrix.identity(r.numRows());
        let p: PermutationMatrix = null;
        if (this.params.enablePivoting)
            p = PermutationMatrix.identity(r.numCols(), PermutationType.Col);

        let rank = 0;
        // make zeroes for each column from the bottom to rowIdx - 1
        for (let j = 0; j != r.numCols(); ++j) {
            // | a
            // | x
            // | b  
            // | 0  /|\
            // | 0 / | \
            if (this.params.enablePivoting) {
                let pivotIdx = j;
                let normSqr = 0.0;
                for (let col = j; col < r.numCols(); ++col) {
                    let curNormSqr = 0;
                    for (let row = j; row < r.numRows(); ++row)
                        curNormSqr += Math.pow(r.get(row, col), 2);
                    if (curNormSqr > normSqr) {
                        pivotIdx = col;
                        normSqr = curNormSqr;
                    }
                }
                if (normSqr < this.params.normSqrTolerance)
                    break;
                if (pivotIdx != j) {
                    r.swapColumns(j, pivotIdx);
                    p.swap(j, pivotIdx);
                }
            }
            for (let i = r.numRows() - 1; i > j; --i) {
                let a = r.get(j, j);
                let b = r.get(i, j);
                let givensCoeffs = givens(a, b);
                // apply to r
                applyGivensFromLeft(r, givensCoeffs, i, j);
                //this.r.set(j, j, givensCoeffs.r);
                //this.r.set(i, j, 0.0);
                // apply to q
                applyTransposeGivensFromRight(q, givensCoeffs, i, j);
            }
            if (!this.params.enablePivoting) {
                if (Math.pow(r.get(j, j), 2) < this.params.normSqrTolerance)
                    return;
            }
            ++rank;
        }
        this.rank = rank;
        this.qDet = 1.0;
        this.q = q;
        this.p = p;
        this.r = r;
    }
    private factorizeHousholder() {
        let r: Matrix = this.r;
        this.r = null;
        let q: Matrix = Matrix.identity(r.numRows());
        let p: PermutationMatrix = null;
        if (this.params.enablePivoting)
            p = PermutationMatrix.identity(r.numCols(), PermutationType.Col);

        let rank = 0;
        // numRows cannot be less than numCols, when numRows == numCols last column should be skipped
        for (let step = 0; step < r.numCols(); ++step) {
            let normSqr = 0.0;
            if (this.params.enablePivoting) {
                let pivotIdx = step;
                for (let col = step; col < r.numCols(); ++col) {
                    let curNormSqr = 0;
                    for (let row = step; row < r.numRows(); ++row)
                        curNormSqr += Math.pow(r.get(row, col), 2);
                    if (curNormSqr > normSqr) {
                        pivotIdx = col;
                        normSqr = curNormSqr;
                    }
                }
                if (normSqr <= this.params.normSqrTolerance)
                    break;
                r.swapColumns(step, pivotIdx);
                p.swap(step, pivotIdx);
            }
            else {
                for (let row = step; row < r.numRows(); ++row)
                    normSqr += Math.pow(r.get(row, step), 2);
                if (normSqr <= this.params.normSqrTolerance)
                    return;
            }
            rank++;
            if (step == r.numRows() - 1) break;
            let v = r.subColumn(step, step, r.numRows() - step);
            let ro = -sign(v.get(0));
            let norm = Math.sqrt(normSqr);
            let firstElement = v.get(0);
            v.set(0, v.get(0) - ro * norm);
            v.scaleSelf(1.0 / Math.sqrt(normSqr - firstElement * firstElement + v.get(0) * v.get(0)));

            // R = (Q3Q2Q1)^A = QT*A
            // QT = Q3(Q2(Q1*I))
            // Q = ((I*Q1)Q2)Q3)
            // QNT = QN because householder matrix is symmetric
            applyHouseholderFromLeft(v, r, step);
            applyHouseholderFromRight(v, q, step);
        }
        // shift signs in order to produce positive diagonal
        for (let step = 0; step < r.numCols(); ++step) {
            for (let k = 0; k < r.numRows(); ++k) {
                q.set(k, step, -q.get(k, step));
                r.set(step, k, -r.get(step, k));
            }
            if (step == r.numRows() - 1) break;
        }
        this.qDet = 1;//(rank & 1 ? -1.0 : 1.0);
        this.rank = rank;
        this.q = q;
        this.p = p;
        this.r = r;
    }
    // Calc A=QR if numRows >= numCols, otherwise AT=QR is computed
    public factorize(A: Matrix | null) {
        this.qDet = 0.0;
        this.q = null;
        this.r = null;
        this.rank = 0;
        this.A = A;
        if (this.A == null) return;

        this.r = this.A.clone();
        if (this.A.numCols() > this.A.numRows())
            this.r.transposeInPlace();
        switch (this.params.zeroingMethod) {
            case ZeroingMethod.Givens:
                this.factorizeGivens();
                break;
            case ZeroingMethod.Housholder:
                this.factorizeHousholder();
                break;
            default:
                assertFail("Invalid value");
        }
        // todo: calculate rank
        if (this.r != null && this.params.makeCompact) {
            this.r.shrinkRows(this.rank);
            this.q.shrinkCols(this.rank);
        }
    }
    constructor(A: Matrix | null, params: OrthogonalDecompositionParams = new OrthogonalDecompositionParams()) {
        this.params = params;
        this.factorize(A);
    }
    get Q(): Matrix {
        return this.q;
    }
    get R(): Matrix {
        return this.r;
    }
    get P(): PermutationMatrix {
        return this.p;
    }
    private solveMatrix(rhs: Matrix): Matrix {
        throw new Error("Not implemented");
    }
    private solveVector(rhs: Vector): Vector {
        assert(this.A.numRows() == rhs.size(), "Incompatible sizes");
        if (this.A.numRows() >= this.A.numCols()) {
            // overdetermined system
            // Ax = b -> QRx = b -> Rx = QT b
            // x = R_1^-1(Q_1^T * b), solve by back substitution
            let b = Matrix.preMulVec(rhs, this.Q);
            let x = Vector.empty(this.A.numCols());
            // back substitution
            for (let row = this.r.numCols() - 1; row >= 0; --row) {
                let sum = b.get(row);
                for (let col = row + 1; col < this.r.numCols(); ++col)
                    sum -= this.R.get(row, col) * x.get(col);
                x.set(row, sum / this.R.get(row, row));
            }
            return x;
        }
        else if (this.A.numCols() > this.A.numRows()) {
            // underdetermined system
            // AT=QR, x = Q(RT^-1 b)
            // forward substitution RTx = b
            // !!! this is not least square norm solution
            let x = Vector.empty(this.A.numCols());
            for (let row = 0; row < this.A.numRows(); ++row) {
                let sum = rhs.get(row);
                for (let col = 1; col < row; ++col)
                    sum -= this.R.get(col, row) * x.get(col);
                x.set(row, sum / this.R.get(row, row));
            }
            x = Matrix.postMulVec(this.Q, x);
            return x;
        }
    }
    solve(rhs: Matrix | Vector): Matrix | Vector {
        if (rhs instanceof Matrix)
            return this.solveMatrix(rhs);
        return this.solveVector(rhs);
    }
    // returns inverse or pseudoinverse depending on the rank of the system
    inverse(): Matrix {
        let inverse: Matrix = Matrix.empty(this.A.numCols(), this.A.numRows());
        if (this.A.numRows() >= this.A.numCols()) {
            for (let row = this.A.numCols() - 1; row >= 0; --row) {
                for (let col = 0; col < this.A.numRows(); ++col) {
                    let sum = this.Q.get(col, row);
                    for (let idx = row + 1; idx < this.A.numCols(); ++idx)
                        sum -= this.R.get(row, idx) * inverse.get(idx, col);
                    inverse.set(row, col, sum / this.R.get(row, row));
                }
            }
            // x = R_1^-1 * Q_1^T 
            // solve R_1 x = Q_1^T
        } else {
            // RT A^-1 = Q
            /*
            for (let row = 0; row < this.A.numRows(); ++row) {
                for (let col = 0; col <= row; ++col) {
                    let rhs = (col == row ? 1 : 0);
                    for (let idx = col; idx < row; ++idx)
                        rhs -= inverse.get(idx, col) * this.R.get(idx, row);
                    inverse.set(row, col, rhs / this.R.get(row, row));
                }
            }
            inverse = Matrix.mul(this.Q, inverse);
            */
            // Alternative
            for (let row = this.A.numRows() - 1; row >= 0; --row) {
                for (let col = 0; col < this.A.numCols(); ++col) {
                    let sum = this.Q.get(col, row);
                    for (let idx = row + 1; idx < this.A.numRows(); ++idx)
                        sum -= this.R.get(row, idx) * inverse.get(col, idx);
                    inverse.set(col, row, sum / this.R.get(row, row));
                }
            }
        }
        return inverse;
    }
    determinant(): any {
        let result = this.qDet;
        for (let i = 0; i < Math.min(this.r.numCols(), this.r.numCols()); ++i)
            result *= this.r.get(i, i);
        return result;
    }
    determinantSign(): number {
        return this.qDet;
    }
}

export class QRTest {
    private A: Matrix = null;
    private q: Matrix = null;
    private r: Matrix = null;
    private rank: number = 0;
    protected p: PermutationMatrix = null;
    private params: OrthogonalDecompositionParams;

    get Params(): OrthogonalDecompositionParams {
        return this.params;
    }
    get Rank(): number {
        return this.rank;
    }
    private factorizeHousholder() {
        let r: Matrix = this.r;
        this.r = null;
        let q: Matrix = Matrix.identity(r.numRows());
        let p: PermutationMatrix = null;
        if (this.params.enablePivoting)
            p = PermutationMatrix.identity(r.numCols(), PermutationType.Col);

        let rank = 0;
        // numRows cannot be less than numCols, when numRows == numCols last column should be skipped
        for (let step = 0; step < r.numCols(); ++step) {
            let normSqr = 0.0;
            if (this.params.enablePivoting) {
                let pivotIdx = step;
                for (let col = step; col < r.numCols(); ++col) {
                    let curNormSqr = 0;
                    for (let row = step; row < r.numRows(); ++row)
                        curNormSqr += Math.pow(r.get(row, col), 2);
                    if (curNormSqr > normSqr) {
                        pivotIdx = col;
                        normSqr = curNormSqr;
                    }
                }
                if (normSqr <= this.params.normSqrTolerance)
                    break;
                r.swapColumns(step, pivotIdx);
                p.swap(step, pivotIdx);
            }
            else {
                for (let row = step; row < r.numRows(); ++row)
                    normSqr += Math.pow(r.get(row, step), 2);
                if (normSqr <= this.params.normSqrTolerance)
                    return;
            }
            rank++;
            if (step == r.numRows() - 1) break;
            let v = r.subColumn(step, step, r.numRows() - step);
            let ro = -sign(v.get(0));
            let norm = Math.sqrt(normSqr);
            let firstElement = v.get(0);
            v.set(0, v.get(0) - ro * norm);
            v.scaleSelf(1.0 / Math.sqrt(normSqr - firstElement * firstElement + v.get(0) * v.get(0)));

            // R = (Q3Q2Q1)^A = QT*A
            // QT = Q3(Q2(Q1*I))
            // Q = ((I*Q1)Q2)Q3)
            // QNT = QN because householder matrix is symmetric
            applyHouseholderFromLeft(v, r, step);
            applyHouseholderFromRight(v, q, step);
        }
        // shift signs in order to produce positive diagonal
        for (let step = 0; step < rank; ++step) {
            for (let k = 0; k < q.numRows(); ++k)
                q.set(k, step, -q.get(k, step));
            for (let k = 0; k < r.numCols(); ++k)
                r.set(step, k, -r.get(step, k));
            if (step == r.numRows() - 1) break;
        }
        this.rank = rank;
        this.q = q;
        this.p = p;
        this.r = r;
    }
    // Calc A=QR if numRows >= numCols, otherwise AT=QR is computed
    public factorize(A: Matrix | null) {
        this.q = null;
        this.r = null;
        this.rank = 0;
        this.A = A;
        if (this.A == null) return;

        this.r = this.A.clone();
        this.factorizeHousholder();
    }
    constructor(A: Matrix | null, params: OrthogonalDecompositionParams = new OrthogonalDecompositionParams()) {
        this.params = params;
        this.factorize(A);
    }
    get Q(): Matrix {
        return this.q;
    }
    get R(): Matrix {
        return this.r;
    }
    get P(): PermutationMatrix {
        return this.p;
    }
    private solveMatrix(rhs: Matrix): Matrix {
        throw new Error("Not implemented");
    }
    private solveVector(rhs: Vector): Vector {
        assert(this.A.numRows() == rhs.size(), "Incompatible sizes");
        let b = Matrix.preMulVec(rhs, this.Q);
        let x = Vector.empty(this.A.numCols());
        for (let row = this.rank - 1; row >= 0; --row) {
            let sum = b.get(row);
            for (let col = row + 1; col < this.r.numCols(); ++col)
                sum -= this.R.get(row, col) * x.get(col);
            x.set(row, sum / this.R.get(row, row));
        }
        return x;
        throw new Error("Not implemented");
        if (this.A.numRows() >= this.A.numCols()) {
            // overdetermined system
            // Ax = b -> QRx = b -> Rx = QT b
            // x = R_1^-1(Q_1^T * b), solve by back substitution
            let b = Matrix.preMulVec(rhs, this.Q);
            let x = Vector.empty(this.A.numCols());
            // back substitution
            for (let row = this.r.numCols() - 1; row >= 0; --row) {
                let sum = b.get(row);
                for (let col = row + 1; col < this.r.numCols(); ++col)
                    sum -= this.R.get(row, col) * x.get(col);
                x.set(row, sum / this.R.get(row, row));
            }
            return x;
        }
        else if (this.A.numCols() > this.A.numRows()) {
            // underdetermined system
            // AT=QR, x = Q(RT^-1 b)
            // forward substitution RTx = b
            // !!! this is not least square norm solution
            let x = Vector.empty(this.A.numCols());
            for (let row = 0; row < this.A.numRows(); ++row) {
                let sum = rhs.get(row);
                for (let col = 1; col < row; ++col)
                    sum -= this.R.get(col, row) * x.get(col);
                x.set(row, sum / this.R.get(row, row));
            }
            x = Matrix.postMulVec(this.Q, x);
            return x;
        }
    }
    solve(rhs: Matrix | Vector): Matrix | Vector {
        if (rhs instanceof Matrix)
            return this.solveMatrix(rhs);
        return this.solveVector(rhs);
    }
    // returns inverse or pseudoinverse depending on the rank of the system
    inverse(): Matrix {
        let inverse: Matrix = Matrix.empty(this.A.numCols(), this.A.numRows());
        for (let row = this.rank - 1; row >= 0; --row) {
            for (let col = 0; col < this.A.numRows(); ++col) {
                let sum = this.Q.get(col, row);
                for (let idx = row + 1; idx < this.A.numCols(); ++idx)
                    sum -= this.R.get(row, idx) * inverse.get(idx, col);
                inverse.set(row, col, sum / this.R.get(row, row));
            }
        }
        // x = R_1^-1 * Q_1^T 
        // solve R_1 x = Q_1^T
        return inverse;
    }
    determinant(): number {
        let result = 1;
        for (let i = 0; i < this.rank; ++i)
            result *= this.r.get(i, i);
        return result;
    }
}
