import { binomial } from "../utils";


test.each([{ n: 5, k: 3, res: 10 }/*, { n: 4, k: 2, res: 6 }, { n: 4, k: 3, res: 4 }, { n: 4, k: 1, res: 4 }, { n: 10, k: 3, res: 120 }*/])('binomial coefficient', (pair: { n: number, k: number, res: number }) => {
    expect(binomial(pair.n, pair.k)).toBeCloseTo(pair.res);
});