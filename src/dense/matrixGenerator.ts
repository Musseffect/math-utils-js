import { complex } from "../complex";
import Matrix from "../denseMatrix";
import { randomNormalDistr } from "../random/generator";
import { QR, ZeroingMethod } from "../solvers/linear systems/qr";
import { randomArray, sign } from "../utils";

export enum SignType {
    Positive,
    Negative,
    Random,
    Ignore
}
// todo: test
export class MatrixGenerator {
    static random(numRows: number, numCols: number): Matrix {
        let data = [];
        for (let i = 0; i < numRows * numCols; i++)
            data.push(Math.random());
        return new Matrix(data, numRows, numCols);
    }
    static randomSymmetric(size: number): Matrix {
        let matrix = Matrix.empty(size, size);
        for (let i = 0; i < size; ++i) {
            for (let j = i; j < size; ++j) {
                let value = Math.random();
                matrix.set(i, j, value);
                matrix.set(j, i, value)
            }
        }
        return matrix;
    }
    static randomOrthogonal(size: number): Matrix {
        let m = Matrix.generate(size, size, (r: number, c: number) => { return randomNormalDistr(); });
        const qrSolver = new QR(m, ZeroingMethod.Housholder, false);
        let result = qrSolver.Q;
        for (let row = 0; row < size; ++row) {
            for (let col = 0; col < size; ++col)
                result.set(row, col, result.get(row, col) * sign(qrSolver.R.get(col, col)));
        }
        return result;
    }
    static randomDiagonal(size: number, shift: number): Matrix {
        let matrix = Matrix.empty(size, size);
        for (let row = 0; row < size; ++row) {
            for (let col = Math.max(0, row - shift); col <= Math.min(size - 1, row + shift); ++col)
                matrix.set(row, col, Math.random());
        }
        return matrix;
    }
    static randomWithEigenvalues(size: number, minEig: number, maxEig: number, makeSymmetric: boolean, sign: SignType): Matrix {
        let eigenvalues = randomArray(size, minEig, maxEig);
        if (sign != SignType.Ignore) {
            for (let i = 0; i < size; ++i) {
                let value = Math.abs(eigenvalues[i]);
                switch (sign) {
                    case SignType.Negative:
                        eigenvalues[i] = -value;
                        break;
                    case SignType.Positive:
                        eigenvalues[i] = value;
                        break;
                    case SignType.Random:
                        eigenvalues[i] = Math.sign(Math.random() - 0.5) * value;
                }
            }
        }
        return this.randomFromEigenvalues(eigenvalues, makeSymmetric);
    }
    static randomFromEigenvalues(eigenvalues: number[], makeSymmetric: boolean): Matrix {
        let M = Matrix.diag(eigenvalues);
        if (!makeSymmetric) {
            for (let row = 0; row + 1 < M.numRows(); ++row) {
                for (let col = row + 1; col < M.numCols(); ++col)
                    M.set(row, col, Math.random() * 2 - 1);
            }
        }
        let Q = this.randomOrthogonal(eigenvalues.length);
        return Matrix.mul(Matrix.mul(Q, M), Q.transpose());
    }
    /**
     * Generate matrix from list of eigenvalues, if eigenvalue is complex it will be treated as conjugate pairs,
     * if it's real it will be counted as one eigenvalue
     * @param eigenvalues 
     */
    static randomFromComplexPairsEigenvalues(eigenvalues: complex[]): Matrix {
        let size = 0;
        for (const eigenvalue of eigenvalues) {
            size++;
            if (Math.abs(eigenvalue.y) != 0)
                size++;
        }
        let m = Matrix.empty(size, size);
        throw new Error("Not implemented");
    }
}