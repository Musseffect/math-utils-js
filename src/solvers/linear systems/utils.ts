import Matrix from "../../denseMatrix"
import { binomial } from "../../utils";


export function hilbertMatrix(size: number): Matrix {
    let m = Matrix.empty(size, size);
    for (let i = 0; i < size; ++i) {
        for (let j = 0; j <= i; ++j) {
            let value = 1 / (i + j + 1);
            m.set(j, i, value);
            m.set(i, j, value)
        }
    }
    return m;
}

export function inverseHilbertMatrix(size: number): Matrix {
    let m = Matrix.empty(size, size);
    for (let i = 0; i < size; ++i) {
        for (let j = 0; j <= i; ++j) {
            let value = (i + j) & 1 ? -1 : 1;
            value *= i + j + 1;
            value *= binomial(size + i, size - j - 1);
            value *= binomial(size + j, size - i - 1);
            value *= Math.pow(binomial(i + j, i), 2);
            m.set(j, i, value);
            m.set(i, j, value)
        }
    }
    return m;
}

/* expect(Matrix.lInfDist(inverseHilbertMatrix(5), new Matrix([25, -300, 1050, -1400, 630, -300, 4800, -18900, 26880, -12600, 1050, -18900, 79380, -117600, 56700,
-1400,, 26880, -117600, 179200, -88200, 630, -12600, 56700, -88200, 44100],5,5)).toBeLessThan(SmallTolerance);
expect(Matrix.lInfDist(hilbertMatrix(5),new Matrix([1, 1/2,1/3,1/4,1/5, 1/2,1/3,1/4,1/5,1/6,1/3,1/4,1/5,1/6,1/7,1/4,1/5,1/6,1/7,1/8,1/5,1/6,1/7,1/8,1/9],5,5))).toBeLessThan(SmallTolerance);
*/