import Matrix from "../denseMatrix";
import { PermutationType, PermutationMatrix } from "../permutationMatrix";
import { SparseMatrixCSR } from "../sparseMatrix";
import { SparseVector } from "../sparseVector";
import Triplet from "../triplet";
import { SmallTolerance, Tolerance, assert } from "../utils";
import Vector from "../vector";

const singularTrivialMatrixTriplets: Triplet[] = [{ row: 1, column: 1, value: 1 },
{ row: 2, column: 2, value: 1 }, { row: 0, column: 0, value: 2 }, { row: 3, column: 3, value: 1 }, { row: 5, column: 5, value: 1 }];
const singularSparseTrivialMatrix = SparseMatrixCSR.fromTriplets(singularTrivialMatrixTriplets, 6, 6);
const singularDenseTrivialMatrix = Matrix.fromTriplets(singularTrivialMatrixTriplets, 6, 6);

const singularNonTrivialMatrixTriplets: Triplet[] = [
    { row: 0, column: 0, value: 1 },
    { row: 0, column: 2, value: 3 },
    { row: 0, column: 4, value: 1 },
    { row: 1, column: 0, value: 2 },
    { row: 1, column: 1, value: 2 },
    { row: 1, column: 5, value: 2 },
    { row: 2, column: 0, value: 2 },
    { row: 2, column: 1, value: 2 },
    { row: 2, column: 2, value: 3 },
    { row: 2, column: 3, value: 1 },
    { row: 2, column: 5, value: 3 },
    { row: 3, column: 4, value: 4 },
    { row: 4, column: 0, value: 3 },
    { row: 4, column: 1, value: -3 },
    { row: 4, column: 2, value: -1 },
    { row: 4, column: 4, value: 2 },
    { row: 5, column: 1, value: 2 },
    { row: 5, column: 3, value: 2 },
    { row: 5, column: 5, value: 4 }
];
const singularSparseNonTrivialMatrix = SparseMatrixCSR.fromTriplets(singularNonTrivialMatrixTriplets, 6, 6);
const singularDenseNonTrivialMatrix = Matrix.fromTriplets(singularNonTrivialMatrixTriplets, 6, 6);

const nonSingularMatrixTriplets: Triplet[] = [{ row: 0, column: 0, value: 1 }, { row: 0, column: 2, value: 3 }, { row: 0, column: 4, value: 1 },
{ row: 1, column: 0, value: 2 }, { row: 1, column: 1, value: 2 }, { row: 1, column: 5, value: 2 }, { row: 2, column: 1, value: 2 }, { row: 2, column: 3, value: 1 }
    , { row: 2, column: 5, value: 3 }, { row: 3, column: 4, value: 4 }, { row: 4, column: 0, value: 3 }, { row: 4, column: 1, value: -3 }, { row: 4, column: 2, value: -1 }
    , { row: 4, column: 4, value: 2 }, { row: 5, column: 1, value: 2 }, { row: 5, column: 3, value: 2 }, { row: 5, column: 5, value: 4 }];
const nonSingularSparseMatrix = SparseMatrixCSR.fromTriplets(nonSingularMatrixTriplets, 6, 6);
const nonSingularDenseMatrix = Matrix.fromTriplets(nonSingularMatrixTriplets, 6, 6);
const nonSingularMatrixDeterminant = 144.0;


test.skip("General sparse matrix", () => {
    expect(Matrix.lInfDistance(singularSparseTrivialMatrix.toDense(), singularDenseTrivialMatrix)).toBeCloseTo(0.0);
    expect(Matrix.lInfDistance(singularSparseNonTrivialMatrix.toDense(), singularDenseNonTrivialMatrix)).toBeCloseTo(0.0);
    expect(Matrix.lInfDistance(nonSingularSparseMatrix.toDense(), nonSingularDenseMatrix)).toBeCloseTo(0.0);

    expect(SparseMatrixCSR.near(SparseMatrixCSR.identity(10).inverse(), SparseMatrixCSR.identity(10))).toBeTruthy();
    let permutationMatrix = new PermutationMatrix([1, 6, 8, 2, 5, 4, 9, 3, 0, 7], PermutationType.Row);
    assert(permutationMatrix.isValid(), "Invalid permutation");
    expect(SparseMatrixCSR.near(permutationMatrix.toSparseMatrix().inverse(), permutationMatrix.inverse().toSparseMatrix()));
    //expect(Math.abs(permutationMatrix.toSparseMatrix().determinant())).toBeCloseTo(1.0);
    //expect(nonSingularSparseMatrix.determinant()).not.toBeCloseTo(0.0);
    expect(SparseMatrixCSR.near(SparseMatrixCSR.mul(nonSingularSparseMatrix.inverse(), nonSingularSparseMatrix), SparseMatrixCSR.identity(nonSingularDenseMatrix.width())));
    // test trivial singular matrix
    //expect(singularSparseTrivialMatrix.determinant()).toBeCloseTo(0.0);
    //expect(singularSparseNonTrivialMatrix.determinant()).toBeCloseTo(0.0);
});

describe('Vector operations', () => {
    test('Construction', () => {
        let dense: Vector = new Vector([0, -1, 2, -1e-8, 0.0, 3, -5, .0]);
        let v: SparseVector = SparseVector.fromVector(dense.data, SmallTolerance);
        expect(Vector.near(v.toDense(), dense, SmallTolerance)).toBeTruthy();
        expect(v.isIndexPresent(0)).toBeFalsy();
        expect(v.isIndexPresent(1)).toBeTruthy();
        expect(v.isIndexPresent(2)).toBeTruthy();
        expect(v.isIndexPresent(3)).toBeFalsy();
        expect(v.isIndexPresent(4)).toBeFalsy();
        expect(v.isIndexPresent(5)).toBeTruthy();
        expect(v.isIndexPresent(6)).toBeTruthy();
        expect(v.isIndexPresent(7)).toBeFalsy();
        for (let i = 0; i < dense.size(); ++i)
            expect(v.get(i)).toBeCloseTo(dense.get(i));

        let v2: SparseVector = SparseVector.empty(8);
        v2.set(2, 2);
        v2.set(1, -1);
        v2.set(5, 3);
        v2.set(6, -5);
        expect(SparseVector.near(v, v2, Tolerance)).toBeTruthy();
        for (let i = 0; i + 1 < v2.elements.length; ++i)
            expect(v2.elements[i].index).toBeLessThan(v2.elements[i + 1].index);

        v.set(1, 0);
        expect(v.isIndexPresent(1)).toBeFalsy();
        v.set(1, -1);
        expect(v.isIndexPresent(1)).toBeTruthy();
        expect(v.get(1)).toBeCloseTo(-1);
        v.set(2, 5);
        expect(v.isIndexPresent(2)).toBeTruthy();
        expect(v.get(2)).toBeCloseTo(5);
        v.set(3, -2);
        expect(v.isIndexPresent(3)).toBeTruthy();
        expect(v.get(3)).toBeCloseTo(-2);

    });
    let d1: Vector = new Vector([0, -1, 2, 0.0, 0.001, 3, -5, .0]);
    let v1: SparseVector = SparseVector.fromVector(d1.data, 0);
    let d2: Vector = new Vector([1, 10.5, 3, 0, -2, 0, 0, 1]);
    let v2: SparseVector = SparseVector.fromVector(d2.data, 0);
    test('Scalar product', () => {
        expect(SparseVector.dot(v1, v2)).toBeCloseTo(Vector.dot(d1, d2));
    });
    test('Pointwise addition', () => {
        expect(Vector.near(SparseVector.add(v1, v2).toDense(), Vector.add(d1, d2))).toBeTruthy();
    });
    test('Pointwise subtraction', () => {
        expect(Vector.near(SparseVector.sub(v1, v2).toDense(), Vector.sub(d1, d2))).toBeTruthy();
    });
    test('Pointwise product', () => {
        expect(Vector.near(SparseVector.mul(v1, v2).toDense(), Vector.mul(d1, d2))).toBeTruthy();
    });
    test('Norms', () => {
        expect(v1.l2Norm()).toBeCloseTo(d1.l2Norm());
        expect(v1.l1Norm()).toBeCloseTo(d1.l1Norm());
        expect(v1.lInfNorm()).toBeCloseTo(d1.lInfNorm());
        expect(v1.squaredLength()).toBeCloseTo(d1.squaredLength());
    });
})

describe('Matrix operations', () => {
    let triplets = [
        { column: 0, row: 0, value: 1 },
        { column: 2, row: 1, value: 1 },
        { column: 4, row: 4, value: 2 },
        { column: 8, row: 1, value: 4 },
        { column: 1, row: 7, value: 7 },
        { column: 1, row: 1, value: 8 },
        { column: 9, row: 7, value: 3 },
        { column: 3, row: 9, value: 8 },
        { column: 3, row: 2, value: 7 }];
    let denseMat = Matrix.fromTriplets(triplets, 10, 10);
    let sparseMat = SparseMatrixCSR.fromTriplets(triplets, 10, 10, 0);
    test('Construction', () => {
        for (let { column, row, value } of triplets)
            expect(sparseMat.get(row, column)).toBeCloseTo(value);
        for (let row = 0; row < 10; ++row) {
            for (let col = 0; col < 10; ++col)
                expect(sparseMat.get(row, col)).toBeCloseTo(denseMat.get(row, col));
        }
        expect(sparseMat.isValid()).toBeTruthy();
        // todo: test order of indices
    });
    //let vec = SparseVector();
    test.skip('Vector multiplication', () => {
    });
    let m2 = SparseMatrixCSR.fromTriplets([], 10, 3, 0);
    test('Multiplication', () => {
        let denseResult = Matrix.mul(sparseMat.toDense(), m2.toDense());
        let sparseResult = SparseMatrixCSR.mul(sparseMat, m2);
        expect(Matrix.lInfDistance(denseResult, sparseResult.toDense())).toBeLessThan(SmallTolerance);
    });
    test.skip('Entrywise product', () => {
        /*let denseResult = Matrix.pointwiseProduct(m1.toDense(), m2.toDense());
        let sparseResult = SparseMatrixCSR.entrywiseProduct(m1, m2);
        expect(Matrix.lInfDistance(denseResult, sparseResult.toDense())).toBeLessThan(SmallTolerance);*/
    });
    test.skip('Kronecker product', () => {

    });
    test('Transpose', () => {
        const sparse = sparseMat;
        let transpose = sparse.transpose();
        expect(transpose.numCols()).toBe(sparse.numRows());
        expect(transpose.numRows()).toBe(sparse.numCols());
        for (let i = 0; i < transpose.numCols(); ++i) {
            for (let j = 0; j < transpose.numRows(); ++j)
                expect(transpose.get(j, i)).toBe(sparse.get(i, j));
        }
        expect(transpose.isValid()).toBeTruthy();
    });
    test.skip('Determinant', () => {

    });
});