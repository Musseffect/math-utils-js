

class Bisection {
    static solve(f: (x: number) => number, a: number, b: number, numIters: number): number {
        let fa = f(a);
        let fb = f(b);
        for (let i = 0; i < numIters; ++i) {
            let x = (a + b) * 0.5;
            let fx = f(x);
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