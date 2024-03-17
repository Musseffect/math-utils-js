import Matrix from "../denseMatrix";
import JSGenerator from "../random/js";
import { SmallTolerance } from "../utils";
import { MatrixGenerator } from "./matrixGenerator";


describe('Random matrix generation tests', () => {
    let generator = new MatrixGenerator(new JSGenerator());
    const MatrixSize = 5;
    test.skip('Random eigenvlaues symmetric', () => {
        // TODO:
    });
    test.skip('Random eigenvlaues general', () => {
        // TODO:
    });
    test('Random orthogonal', () => {
        let matrix = generator.randomOrthogonal(MatrixSize);
        expect(matrix.isOrthogonal()).toBeTruthy();
    });
    test('Random symmetric', () => {
        let matrix = generator.randomSymmetric(MatrixSize);
        expect(matrix.isSymmetric()).toBeTruthy();
    });
});

describe('Special matrix generators', () => {
    test('Hilbert matrix', () => {
        const hilbert = new Matrix(
            [
                1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
                1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6,
                1 / 3, 1 / 4, 1 / 5, 1 / 6, 1 / 7,
                1 / 4, 1 / 5, 1 / 6, 1 / 7, 1 / 8,
                1 / 5, 1 / 6, 1 / 7, 1 / 8, 1 / 9], 5, 5);
        const invHilbert = new Matrix(
            [
                25, -300, 1050, -1400, 630,
                -300, 4800, -18900, 26880, -12600,
                1050, -18900, 79380, -117600, 56700,
                -1400, 26880, -117600, 179200, -88200,
                630, -12600, 56700, -88200, 44100], 5, 5);
        expect(Matrix.lInfDistance(MatrixGenerator.inverseHilbertMatrix(5), invHilbert)).toBeLessThan(SmallTolerance);
        expect(Matrix.lInfDistance(MatrixGenerator.hilbertMatrix(5), hilbert)).toBeLessThan(SmallTolerance);
    });
});