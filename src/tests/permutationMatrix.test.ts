import Matrix from "../denseMatrix";
import PermutationMatrix from "../permutationMatrix";
import { assert } from "../utils";
import Vector from "../vector";


test("Permutation matrix", () => {
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
    //console.log(`identity: ${Matrix.identity(10).toString()}`);
    //console.log(`rowPermutations: ${rowPermutations.toMatrix().toString()}`);
    const rowPermuted = Matrix.mul(rowPermutations.toMatrix(), Matrix.identity(10));
    //console.log(`rowPermuted: ${rowPermuted.toString()}`);
    const colPermuted = Matrix.mul(Matrix.identity(10), colPermutations.toMatrix());
    const fullPermuted = Matrix.mul(rowPermuted, colPermutations.toMatrix());
    for (let i = 0; i < 10; ++i) {
        const idx1 = rowPermutations.permuteIndex(i, i);
        expect(rowPermuted.get(idx1.row, idx1.column)).toBeCloseTo(1.0);
        expect(rowPermutations.value(idx1.row)).toEqual(i);
        const idx2 = colPermutations.permuteIndex(i, i);
        expect(colPermuted.get(idx2.row, idx2.column)).toBeCloseTo(1.0);
        expect(colPermutations.value(idx2.column)).toEqual(i);
        const idx3 = colPermutations.permuteIndex(i, i);
        const idx4 = rowPermutations.permuteIndex(idx3.row, idx3.column);
        expect(idx4.column).toEqual(idx3.column);
        expect(idx3.row).toEqual(i);
        expect(fullPermuted.get(idx4.row, idx4.column)).toBeCloseTo(1.0);
        expect(colPermutations.value(idx4.column)).toEqual(i);
        expect(rowPermutations.value(idx3.column)).toEqual(i);
    }
    assert(rowPermutations.isValid(), "Invalid permutation");
    assert(colPermutations.isValid(), "Invalid permutation");
    expect(Matrix.near(rowPermuted, rowPermutations.toMatrix())).toBeTruthy();
    expect(Matrix.near(colPermuted, colPermutations.toMatrix())).toBeTruthy();
    expect(Matrix.near(fullPermuted, Matrix.mul(rowPermutations.toMatrix(), colPermutations.toMatrix()))).toBeTruthy();
    expect(Matrix.near(Matrix.mul(colPermutations.inverse().toMatrix(), colPermutations.toMatrix()), Matrix.identity(10))).toBeTruthy();
    expect(Matrix.near(Matrix.mul(rowPermutations.inverse().toMatrix(), rowPermutations.toMatrix()), Matrix.identity(10))).toBeTruthy();

    //expect(Matrix.near(permutationMatrix.toMatrix().inverse(), permutationMatrix.inverse().toMatrix()));
    //expect(Math.abs(permutationMatrix.toMatrix().determinant())).toBeCloseTo(1.0);

    const generatedMat = Matrix.generate(10, 10, (r: number, c: number) => r * 10 + c);
    const generatedVec = Vector.generate(10, (i: number) => i);
    const rowPermutedMat = rowPermutations.permuteMatrix(generatedMat);
    const colPermutedMat = colPermutations.permuteMatrix(generatedMat);
    const fullPermutedMat = colPermutations.permuteMatrix(rowPermutedMat);
    const rowPermutedVec = rowPermutations.permuteVector(generatedVec);
    const colPermutedVec = colPermutations.permuteVector(generatedVec);
    const fullPermutedVec = colPermutations.permuteVector(rowPermutedVec);
    // console.log(`Row permutation: ${rowPermutations.array()}`);
    // console.log(`rowPermutedMat : ${rowPermutedMat.toString()}`);
    // console.log(`rowPermutationsMatrix : ${rowPermutations.toMatrix().toString()}`);
    // console.log(`generatedMatMatrix : ${generatedMat.toString()}`);
    // console.log(`generatedMatMatrix : ${Matrix.mul(rowPermutations.toMatrix(), generatedMat).toString()}`);
    expect(Matrix.lInfDistance(rowPermutedMat, Matrix.mul(rowPermutations.toMatrix(), generatedMat))).toBeCloseTo(0);
    expect(Matrix.lInfDistance(colPermutedMat, Matrix.mul(generatedMat, colPermutations.toMatrix()))).toBeCloseTo(0);
    expect(Matrix.lInfDistance(fullPermutedMat, Matrix.mul(Matrix.mul(rowPermutations.toMatrix(), generatedMat), colPermutations.toMatrix()))).toBeCloseTo(0);
    expect(Vector.lInfDistance(rowPermutedVec, Matrix.postMulVec(rowPermutations.toMatrix(), generatedVec))).toBeCloseTo(0);
    expect(Vector.lInfDistance(colPermutedVec, Matrix.preMulVec(generatedVec, colPermutations.toMatrix()))).toBeCloseTo(0);
    expect(Vector.lInfDistance(fullPermutedVec, Matrix.preMulVec(Matrix.postMulVec(rowPermutations.toMatrix(), generatedVec), colPermutations.toMatrix()))).toBeCloseTo(0);
});