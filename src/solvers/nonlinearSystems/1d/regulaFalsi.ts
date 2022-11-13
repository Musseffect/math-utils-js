import { assert, SmallEpsilon } from "../../../utils";

class RegulaFalsi {
    static solve(f: (x: number) => number, a: number, b: number, numIters: number, tolerance: number = SmallEpsilon): number {
        let fa = f(a);
        if (Math.abs(fa) < tolerance) return a;
        let fb = f(b);
        if (Math.abs(fb) < tolerance) return b;
        assert(fb * fa <= 0, "Bad initial conditions");
        for (let i = 0; i < numIters; ++i) {
            let x = a - (b - a) / (fb - fa);
            let fx = f(x);
            if (Math.abs(fx) < tolerance)
                return x;
            if (fx * fa < fx * fb) {
                b = x;
                fb = fx;
            } else {
                a = x;
                fa = fx;
            }
        }
        return Math.abs(fa) < Math.abs(fb) ? a : b;
    }
}