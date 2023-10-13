import { complex } from "./complex";
import Matrix from "./denseMatrix";
import { SmallTolerance, SmallestTolerance, assert, assertFail } from "./utils";

// coefficients from 0 to N -> c_0 + c_1 * x + c_2 * x^2 + ... + c_N * x^N
export class Polynomial {
    coeffs: number[];
    constructor(coeffs: number[]) {
        assert(coeffs.length >= 1, "Invalid size");
        this.coeffs = coeffs;
    }
    numCoeffs(): number {
        return this.coeffs.length;
    }
    degree(): number {
        return this.coeffs.length - 1;
    }
    eval(x: number): number {
        let value = 0.0;
        for (let i = this.numCoeffs() - 1; i > 0; --i) {
            value += this.coeffs[i];
            value *= x;
        }
        return value + this.coeffs[0];
    }
    derivative(): Polynomial {
        assert(this.coeffs.length > 1, "Invalid size");
        let result: number[] = [];
        for (let i = 1; i < this.numCoeffs(); ++i) {
            result.push(this.coeffs[i] * i);
        }
        return new Polynomial(result);
    }
    shrink(): void {
        while (this.coeffs.length > 0 && this.coeffs[-1] == 0)
            this.coeffs.pop();
    }
    companionMatrix(): Matrix {
        let matrix = Matrix.empty(this.degree(), this.degree());
        let lastCoeff = this.coeffs[-1];
        for (let i = 0; i < this.degree(); ++i) {
            matrix.set(i + 1, i, 1);
            matrix.set(i, this.degree(), -this.coeffs[i] / lastCoeff)
        }
        return matrix;
    }
}

export class PolynomialSolver {
    protected numIters: number = 10;
    protected tol: number = SmallestTolerance;
    protected delta: number = 1.0;
    constructor(numIters: number = 10, tol: number = SmallestTolerance, delta: number = 1) {
        this.numIters = numIters;
        this.tol = tol;
        this.delta = delta;
    }
    // bisection + newton
    // todo: try false position and fancy false position from matlab paper
    protected findRoot(f: Polynomial, df: Polynomial, xMin: number, xMax: number, fMin: number, fMax: number): number {
        assert(Math.sign(fMin) != Math.sign(fMax), "Incorrect interval");
        // handle case with both infinity bounds, it should only arise when polynomial is odd
        if (xMin == Number.NEGATIVE_INFINITY && xMax == Number.POSITIVE_INFINITY) {
            let x0 = 0.0;
            let value = f.eval(x0);
            if (Math.sign(value) == Math.sign(fMin)) {
                xMin = x0;
                fMin = value;
            } else {
                xMax = x0;
                fMax = value;
            }
        }
        for (let i = 0; i < this.numIters; ++i) {
            // handle infinity on one of the sides
            let xCenter: number;
            if (xMin == Number.NEGATIVE_INFINITY)
                xCenter = xMax - this.delta;
            else if (xMax == Number.POSITIVE_INFINITY)
                xCenter = xMin + this.delta;
            else
                xCenter = (xMax + xMin) * 0.5;

            let value = f.eval(xCenter);
            let derivative = df.eval(xCenter);
            if (Math.abs(derivative) > SmallestTolerance) {
                let xNext = xCenter - value / derivative;
                if (xNext > xMin && xNext < xMax) {
                    xCenter = xNext;
                    value = f.eval(xNext);
                    if (Math.abs(xMin - xNext) <= this.tol) {
                        if (Math.sign(fMin) == Math.sign(value))
                            xCenter += this.tol;
                        else xCenter -= this.tol;
                    }
                }
            }
            if (Math.min(Math.abs(xCenter - xMax), Math.abs(xCenter - xMin)) < this.tol)
                return xCenter;
            if (Math.sign(value) == Math.sign(fMin)) {
                xMin = xCenter;
                fMin = value;
            }
            else { xMax = xCenter; fMax = value; }
        }
        return (xMax + xMin) * 0.5;
    }
    private infSign(x: number, f: Polynomial) {
        for (let i = f.numCoeffs() - 1; i >= 0; --i) {
            if (f.coeffs[i] == 0) continue;
            let sign = Math.sign(f.coeffs[i]);
            if (i & 1)
                return Math.sign(x) * sign;
            return sign;
        }
        assertFail("Empty polynomial");
    }
    // ax + b
    static solveLinear(a: number, b: number): number {
        return -b / a;
    }
    // ax^2 + bx + c
    static solveQuadratic(a: number, b: number, c: number, xStart: number = Number.NEGATIVE_INFINITY, xEnd: number = Number.POSITIVE_INFINITY): number[] {
        let d = b * b - 4 * a * c;
        if (d < 0.0)
            return [];
        d = (b + Math.sign(b) * Math.sqrt(d)) / 2;
        let x0 = -c / d;
        let x1 = -d / a;
        let result: number[] = [];
        if (x0 >= xStart && x0 <= xEnd)
            result.push(x0);
        if (x1 >= xStart && x1 <= xEnd)
            result.push(x1);
        return result
    }
    solveCubic(a: number, b: number, c: number, d: number, deflation: boolean = true, xStart: number = Number.NEGATIVE_INFINITY, xEnd: number = Number.POSITIVE_INFINITY): number[] {
        let f = new Polynomial([d, c, b, a]);
        let df = new Polynomial([c, 2 * b, 3 * a]);
        let criticalPoints = PolynomialSolver.solveQuadratic(df.coeffs[0], df.coeffs[1], df.coeffs[2], xStart, xEnd);
        criticalPoints.unshift(xStart);
        criticalPoints.push(xEnd);
        let roots: number[] = [];
        let i = 0;
        for (; i < criticalPoints.length - 1; ++i) {
            let xMin = criticalPoints[i];
            let xMax = criticalPoints[i + 1];
            let fMin = Number.isFinite(xMin) ? f.eval(xMin) : this.infSign(xMin, f);
            let fMax = Number.isFinite(xMax) ? f.eval(xMax) : this.infSign(xMax, f);
            if (Math.sign(fMin) != Math.sign(fMax)) {
                let root = this.findRoot(f, df, xMin, xMax, fMin, fMax);
                roots.push(root);
                if (!deflation) continue;
                // deflation
                let aq = a;
                let bq = b + aq * root;
                let cq = d;
                let otherRoots = PolynomialSolver.solveQuadratic(aq, bq, cq, xMax, xEnd);
                roots = roots.concat(otherRoots);
                break;
            }
        }
        if (deflation && roots.length > 1) {
            for (let j = i + 1; j < criticalPoints.length; ++j) {
                let xMin = criticalPoints[i];
                let xMax = criticalPoints[i + 1];
                let fMin = Number.isFinite(xMin) ? f.eval(xMin) : this.infSign(xMin, f);
                let fMax = Number.isFinite(xMax) ? f.eval(xMax) : this.infSign(xMax, f);
                if (Math.sign(fMin) != Math.sign(fMax)) {
                    let rootIndex = roots.findIndex((value) => {
                        return value >= xMin && value <= xMax;
                    });
                    if (rootIndex >= 0) {
                        let root = roots[rootIndex];
                        let fRoot = f.eval(root);
                        if (Math.sign(fMin) == Math.sign(fRoot)) {
                            xMin = root;
                            fMin = fRoot;
                        }
                        else {
                            xMax = root;
                            fMax = fRoot;
                        }
                        // refine root
                        root = this.findRoot(f, df, xMin, xMax, fMin, fMax);
                        roots[rootIndex] = root;
                    }
                }
            }
        }
        return roots;
    }
    // a[0] + a[1]*x + ...a[n]*x^n
    static solve(a: number[]): number[] {
        assert(a.length > 1, "Empty coefficients");
        if (a.length == 2)
            return [this.solveLinear(a[1], a[0])];
        if (a.length == 3) return this.solveQuadratic(a[2], a[1], a[0]);
    }
    static eval(x: number, a: number[]) {
        let value = 0.0;
        for (let l = a.length - 1; l > 0; --l) {
            value += a[l];
            value *= x;
        }
        return value + a[0];
    }
    solveInRegion(polynomial: Polynomial, xStart: number = Number.NEGATIVE_INFINITY, xEnd: number = Number.POSITIVE_INFINITY) {
        assert(xStart < xEnd, "Invalid interval");
        polynomial.shrink();
        let stack: Polynomial[] = [polynomial];
        for (let i = polynomial.numCoeffs(); i >= 2; --i)
            stack.push(stack[-1].derivative());

        const quadraticPolynomial = stack[-1];
        assert(quadraticPolynomial.numCoeffs() == 2, "Quadratic polynomial expected");
        // todo: use cubic solver with deflation here
        let criticalPoints: number[] = PolynomialSolver.solveQuadratic(quadraticPolynomial.coeffs[2], quadraticPolynomial.coeffs[1], quadraticPolynomial.coeffs[0], xStart, xEnd);
        criticalPoints.unshift(xStart);
        criticalPoints.push(xEnd);
        while (stack.length != 1) {
            let df = stack[-1];
            let f = stack[-2];
            let newRoots: number[] = [];
            newRoots.push(xStart);
            for (let i = 0; i < criticalPoints.length - 1; ++i) {
                let xMin = criticalPoints[i];
                let xMax = criticalPoints[i + 1];
                let fMin = Number.isFinite(xMin) ? f.eval(xMin) : this.infSign(xMin, f);
                let fMax = Number.isFinite(xMax) ? f.eval(xMax) : this.infSign(xMax, f);
                if (Math.sign(fMin) != Math.sign(fMax)) {
                    let newRoot = this.findRoot(f, df, xMin, xMax, fMin, fMax);
                    newRoots.push(newRoot);
                }
            }
            newRoots.push(xEnd);
            criticalPoints = newRoots;
            stack.pop();
        }
        return criticalPoints;
    }
}
export function generatePolynomialWithComplexRoots(roots: complex[]) {
    assert(roots.length != 0, "Zero roots");
    assert(roots.length <= 63, "Too many roots");
    let coeffs: complex[] = [];
    for (let i = 0; i < roots.length; ++i)
        coeffs.push(complex.empty());
    for (let i = 0; i < 2 << roots.length; ++i) {
        let coeffID = 0;
        let value = new complex(-1.0, 0.0);
        for (let j = 0, k = i; j < roots.length; ++j, k = k >> 1) {
            coeffID += k & 1;
            value = complex.mul(value, k & 1 ? new complex(-1, 0) : roots[j]);
        }
        coeffs[coeffID].addSelf(value);
    }
    let realCoeffs: number[] = [];
    for (const coeff of coeffs) {
        assert(Math.abs(coeff.y) < SmallTolerance, "Complex coefficient");
        realCoeffs.push(coeff.x);
    }
    return new Polynomial(realCoeffs);
}

export function generatePolynomial(roots: number[]): Polynomial {
    assert(roots.length != 0, "Zero roots");
    assert(roots.length <= 63, "Too many roots");
    let coeffs = Array(roots.length);
    for (let i = 0; i < 2 << roots.length; ++i) {
        let coeffID = 0;
        let value = 1.0;
        for (let j = 0, k = i; j < roots.length; ++j, k = k >> 1) {
            coeffID += k & 1;
            value *= k & 1 ? -1 : roots[j];
        }
        coeffs[coeffID] += value;
    }
    return new Polynomial(coeffs);
}