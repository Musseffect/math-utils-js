import Matrix from "../../denseMatrix";
import Vector from "../../vector";
// regular linear
export class PiecewiseLinear {
    static run(x: number[], y: number[], intervals: number, isPeriodic: boolean) {
        let A: Matrix = Matrix.empty(intervals + 1, intervals + 1);
        let B: Vector = Vector.empty(intervals + 1);
        let bins: number[][] = [];
        for (let i = 0; i < intervals; i++)
            bins.push([]);
        let max = x[0];
        let min = x[0];
        for (let i = 0; i < x.length; i++) {
            max = Math.max(max, x[i]);
            min = Math.min(min, x[i])
        }
        let dx = (max - min) / intervals;
        for (let i = 0; i < x.length; i++) {
            let bin = Math.min(Math.floor((x[i] - min) / dx), intervals - 1);
            bins[bin].push(i);
        }
        let a: number[] = [];
        for (let i = 0; i < intervals + 1; i++) {
            let ai = min + dx * i;
            a.push(ai);
            if (i - 1 >= 0) {
                let aPrev = ai - dx;
                let bin = bins[i - 1];
                for (let j = 0; j < bin.length; j++) {
                    let index = bin[j];
                    let _x = x[index];
                    let _y = y[index];
                    B.set(B.get(i) + _y * (_x - aPrev) / dx, i);
                    A.set(A.get(i, i - 1) + (ai - _x) * (_x - aPrev) / (dx * dx), i, i - 1);
                    A.set(A.get(i, i) + (_x - aPrev) * (_x - aPrev) / (dx * dx), i, i);
                }
            }
            if (i + 1 <= intervals) {
                let aNext = ai + dx;
                let bin = bins[i];
                for (let j = 0; j < bin.length; j++) {
                    let index = bin[j];
                    let _x = x[index];
                    let _y = y[index];
                    B.set(B.get(i) + _y * (aNext - _x) / dx, i);
                    A.set(A.get(i, i) + (aNext - _x) * (aNext - _x) / (dx * dx), i, i);
                    A.set(A.get(i, i + 1) + (_x - ai) * (aNext - _x) / (dx * dx), i, i + 1);
                }
            }
        }
        if (isPeriodic) {
            B.set(B.get(0) + B.get(intervals), 0);
            A.set(A.get(0, 0) + A.get(intervals, intervals), 0, 0);
            A.set(A.get(intervals, intervals - 1) + A.get(0, intervals - 1), 0, intervals - 1);
            A.set(A.get(intervals - 1, 0) + A.get(intervals - 1, intervals), intervals - 1, 0);

            A.set(1, intervals, 0);
            A.set(-1, intervals, intervals);
            A.set(0, 0, intervals);
            for (let i = 1; i < intervals; i++) {
                A.set(0, intervals, i);
                A.set(0, i, intervals);
            }
            B.set(0, intervals);
        }
        let b = Matrix.solve(A, B);
        let error = 0;
        for (let i = 0; i < intervals; i++) {
            let ai = min + dx * i;
            let aNext = ai + dx;
            let bin = bins[i];
            for (let j = 0; j < bin.length; j++) {
                let index = bin[j];
                let _x = x[index];
                let _y = y[index];
                error += Math.pow((b.get(i + 1) * (_x - ai) + b.get(i) * (aNext - _x)) / dx - _y, 2);
            }
        }
        return { a: a, b: b.data, error: error };
    }
}