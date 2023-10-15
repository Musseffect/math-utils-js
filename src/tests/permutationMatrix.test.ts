import Matrix from "../denseMatrix";
import PermutationMatrix from "../permutationMatrix";
import { assert } from "../utils";


test.skip("Permutation matrix", () => {
    let rowPermutations = new PermutationMatrix([1, 6, 8, 2, 5, 4, 9, 3, 0, 7], true);
    let colPermutations = new PermutationMatrix([1, 6, 8, 2, 5, 4, 9, 3, 0, 7], false);
    //expect(rowPermutations.determinant()).toBeCloseTo(rowPermutations.toMatrix().determinant());
    //expect(colPermutations.determinant()).toBeCloseTo(colPermutations.toMatrix().determinant());
    /*const mat = Matrix.empty(rowPermutations.size(), rowPermutations.size());
    for (let i = 0; i < mat.width(); ++i)
    {
        for (let j = 0; j < mat.height(); ++j)
            mat.set(j, i, j * mat.width() + i);
    }*/
    const rowPermuted = Matrix.mul(rowPermutations.toMatrix(), Matrix.identity(10));
    const colPermuted = Matrix.mul(Matrix.identity(10), colPermutations.toMatrix());
    const fullPermuted = Matrix.mul(rowPermuted, colPermutations.toMatrix());
    for (let i = 0; i < 10; ++i) {
        const idx1 = rowPermutations.permuteIndex(i, i);
        expect(rowPermuted.get(idx1.row, idx1.column)).toBeCloseTo(1.0);
        const idx2 = colPermutations.permuteIndex(i, i);
        expect(colPermuted.get(idx2.row, idx2.column)).toBeCloseTo(1.0);
        const idx3 = colPermutations.permuteIndex(idx1.row, idx1.column);
        expect(fullPermuted.get(idx3.row, idx3.column)).toBeCloseTo(1.0);
    }
    assert(rowPermutations.isValid(), "Invalid permutation");
    assert(colPermutations.isValid(), "Invalid permutation");
    expect(Matrix.near(rowPermuted, rowPermutations.toMatrix())).toBeTruthy();
    expect(Matrix.near(colPermuted, colPermutations.toMatrix())).toBeTruthy();
    expect(Matrix.near(fullPermuted, Matrix.mul(rowPermutations.toMatrix(), colPermutations.toMatrix()))).toBeTruthy();
    expect(Matrix.near(Matrix.mul(colPermutations.inverse().toMatrix(), colPermutations.toMatrix()), Matrix.identity(10))).toBeTruthy();
    expect(Matrix.near(Matrix.mul(rowPermutations.inverse().toMatrix(), colPermutations.toMatrix()), Matrix.identity(10))).toBeTruthy();

    //expect(Matrix.near(permutationMatrix.toMatrix().inverse(), permutationMatrix.inverse().toMatrix()));
    //expect(Math.abs(permutationMatrix.toMatrix().determinant())).toBeCloseTo(1.0);

});