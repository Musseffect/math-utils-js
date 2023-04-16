
/**
 * 1-dimensional continuous optimization
 */
export abstract class Problem1D {
    min: number;
    max: number;
    initial: number;
    getBounds(): { min: number, max: number } {
        return { min: this.min, max: this.max };
    }
    getInitial(): number {
        return this.initial;
    }
    abstract f(x: number): number;
    abstract dfdx(x: number): number;
    abstract dfdxdx(x: number): number;
}

export interface Solver1D {
    solve(problem: Problem1D): number;
}

export class GoldenSearch1D implements Solver1D {
    epsilon: number;
    iterations: number;
    constructor(epsilon: number, iterations: number) {
        this.epsilon = epsilon;
        this.iterations = iterations;
    }
    solve(problem: Problem1D): number {
        let { min, max } = problem.getBounds();
        let a = min, b = max;
        let t = (1.0 + Math.sqrt(5.0)) / 2.0;

        let x2: number = a + (b - a) / t;
        let x1: number = a + b - x2;
        let fx1: number = problem.f(x2);
        let fx2: number = problem.f(x1);
        let i = 0;
        while (i < this.iterations) {
            if (fx1 < fx2) {
                b = x2;
                x2 = x1;
                fx2 = fx1;
                x1 = a + b - x2;
                fx1 = problem.f(x1);
            }
            else {
                a = x1;
                x1 = x2;
                fx1 = fx2;
                x2 = a + (b - a) / t;
                fx2 = problem.f(x2);
            }
            if (Math.abs(a - b) < this.epsilon) {
                break;
            }
            i++;
        }
        let x = (a + b) * 0.5;
        return x;
    }
}