import { SmallestEpsilon, assert, assertFail } from "../../utils";

export class Polynomial {
    coeffs: number[];
    constructor(coeffs: number[]) {
        assert(coeffs.length > 1, "Invalid size");
    }
    degree(): number {
        return this.coeffs.length;
    }
    eval(x: number): number {
        let value = 0.0;
        for (let i = this.degree() - 1; i > 0; --i) {
            value += this.coeffs[i];
            value *= x;
        }
        return value + this.coeffs[0];
    }
    derivative(): Polynomial {
        let result: number[] = [];
        for (let i = 1; i < this.degree(); ++i) {
            result.push(this.coeffs[i] * i);
        }
        return new Polynomial(result);
    }
    shrink() {
        while (this.coeffs.length > 0 && this.coeffs[-1] == 0)
            this.coeffs.pop();
    }
}

export class PolynomialSolver {
    protected numIters: number = 10;
    protected tol: number = SmallestEpsilon;
    protected delta: number = 1.0;
    constructor(numIters: number = 10, tol: number = SmallestEpsilon, delta: number = 1) {
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
            if (Math.abs(derivative) > SmallestEpsilon) {
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
        for (let i = f.degree() - 1; i >= 0; --i) {
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
    static solveQuadratic(a: number, b: number, c: number): number[] {
        let d = b * b - 4 * a * c;
        if (d < 0.0)
            return [];
        d = (b + Math.sign(b) * Math.sqrt(d)) / 2;
        return [-c / d, -d / a]
    }
    static solveCubic(a: number, b: number, c: number, d: number) {
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
    solveInRegion(polynomial: Polynomial, x0: number = Number.NEGATIVE_INFINITY, x1: number = Number.POSITIVE_INFINITY) {
        assert(x0 < x1, "Invalid interval");
        polynomial.shrink();
        let stack: Polynomial[] = [polynomial];
        for (let i = polynomial.degree(); i >= 2; --i)
            stack.push(stack[-1].derivative());

        const quadraticPolynomial = stack[-1];
        assert(quadraticPolynomial.degree() == 2, "Quadratic polynomial expected");
        // todo: implement cubic solver with deflation and use it here
        let criticalPoints: number[] = PolynomialSolver.solveQuadratic(quadraticPolynomial.coeffs[2], quadraticPolynomial.coeffs[1], quadraticPolynomial.coeffs[0]);
        // remove roots outside of interval
        while (criticalPoints.length > 0 && criticalPoints[-1] > x1)
            criticalPoints.pop();
        while (criticalPoints.length > 0 && criticalPoints[0] < x0)
            criticalPoints.shift();
        criticalPoints.unshift(x0);
        criticalPoints.push(x1);
        while (stack.length != 1) {
            let df = stack[-1];
            let f = stack[-2];
            let newRoots: number[] = [];
            newRoots.push(x0);
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
            newRoots.push(x1);
            criticalPoints = newRoots;
            stack.pop();
        }
        return criticalPoints;
    }
}
