import Matrix from "../denseMatrix";
import PermutationMatrix from "../permutationMatrix";
import SparseMatrix from "../sparseMatrix";
import SparseVector from "../sparseVector";
import Triplet from "../triplet";
import { SmallTolerance, assert } from "../utils";
import Vector from "../vector";

const singularTrivialMatrixTriplets: Triplet[] = [{ row: 1, column: 1, value: 1 },
{ row: 2, column: 2, value: 1 }, { row: 0, column: 0, value: 2 }, { row: 3, column: 3, value: 1 }, { row: 5, column: 5, value: 1 }];
const singularSparseTrivialMatrix = SparseMatrix.fromTriplets(singularTrivialMatrixTriplets, 6, 6);
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
const singularSparseNonTrivialMatrix = SparseMatrix.fromTriplets(singularNonTrivialMatrixTriplets, 6, 6);
const singularDenseNonTrivialMatrix = Matrix.fromTriplets(singularNonTrivialMatrixTriplets, 6, 6);

const nonSingularMatrixTriplets: Triplet[] = [{ row: 0, column: 0, value: 1 }, { row: 0, column: 2, value: 3 }, { row: 0, column: 4, value: 1 },
{ row: 1, column: 0, value: 2 }, { row: 1, column: 1, value: 2 }, { row: 1, column: 5, value: 2 }, { row: 2, column: 1, value: 2 }, { row: 2, column: 3, value: 1 }
    , { row: 2, column: 5, value: 3 }, { row: 3, column: 4, value: 4 }, { row: 4, column: 0, value: 3 }, { row: 4, column: 1, value: -3 }, { row: 4, column: 2, value: -1 }
    , { row: 4, column: 4, value: 2 }, { row: 5, column: 1, value: 2 }, { row: 5, column: 3, value: 2 }, { row: 5, column: 5, value: 4 }];
const nonSingularSparseMatrix = SparseMatrix.fromTriplets(nonSingularMatrixTriplets, 6, 6);
const nonSingularDenseMatrix = Matrix.fromTriplets(nonSingularMatrixTriplets, 6, 6);
const nonSingularMatrixDeterminant = 144.0;


test.skip("Sparse vector", () => {
    let dense: Vector = new Vector([0, -1, 2, -1e-8, 0.0, 3, -5, .0]);
    let v: SparseVector = SparseVector.fromVector(dense.data, SmallTolerance);
    expect(Vector.near(v.toDense(), dense, SmallTolerance)).toBeTruthy();
    expect(v.isNonZero(0)).toBeFalsy();
    expect(v.isNonZero(1)).toBeTruthy();
    expect(v.isNonZero(2)).toBeTruthy();
    expect(v.isNonZero(3)).toBeFalsy();
    expect(v.isNonZero(4)).toBeFalsy();
    expect(v.isNonZero(5)).toBeTruthy();
    expect(v.isNonZero(6)).toBeTruthy();
    expect(v.isNonZero(7)).toBeFalsy();
    for (let i = 0; i < dense.size(); ++i)
        expect(v.get(i)).toBeCloseTo(dense.get(i));

    expect(v.squaredLength()).toBeCloseTo(dense.squaredLength());
    expect(v.l1Norm()).toBeCloseTo(dense.l1Norm());
    expect(v.l2Norm()).toBeCloseTo(dense.l2Norm());
    expect(v.lInfNorm()).toBeCloseTo(dense.lInfNorm());

    v.set(1, 0);
    expect(v.isNonZero(1)).toBeFalsy();
    v.set(1, -1);
    expect(v.isNonZero(1)).toBeTruthy();
    expect(v.get(1)).toBeCloseTo(-1);
    v.set(2, 5);
    expect(v.isNonZero(2)).toBeTruthy();
    expect(v.get(2)).toBeCloseTo(5);
    v.set(3, -2);
    expect(v.isNonZero(3)).toBeTruthy();
    expect(v.get(3)).toBeCloseTo(-2);

    v = SparseVector.fromVector(dense.data, SmallTolerance);
    let dense2: Vector = new Vector([3, 0.1, 0.0, 1e-7, 0.3, -2.0, -6.0, 7.0]);
    let v2 = SparseVector.fromVector(dense2.data, SmallTolerance);
    expect(Vector.near(v2.toDense(), dense2, SmallTolerance));
    console.log(SparseVector.add(v, v2).toDense().toString());
    console.log(Vector.add(dense, dense2).toString());
    expect(Vector.near(SparseVector.add(v, v2).toDense(), Vector.add(dense, dense2))).toBeTruthy();
    expect(Vector.near(SparseVector.sub(v, v2).toDense(), Vector.sub(dense, dense2))).toBeTruthy();
    expect(Vector.near(SparseVector.mul(v, v2).toDense(), Vector.mul(dense, dense2))).toBeTruthy();
    expect(SparseVector.dot(v, v2),).toBeCloseTo(Vector.dot(dense, dense2));

    //expect(Vector.near(SparseVector.div(v, v2).toDense(), Vector.div(dense, dense2))).toBeTruthy();
});

test.skip("General sparse matrix", () => {



    expect(Matrix.lInfDistance(singularSparseTrivialMatrix.toDense(), singularDenseTrivialMatrix)).toBeCloseTo(0.0);
    expect(Matrix.lInfDistance(singularSparseNonTrivialMatrix.toDense(), singularDenseNonTrivialMatrix)).toBeCloseTo(0.0);
    expect(Matrix.lInfDistance(nonSingularSparseMatrix.toDense(), nonSingularDenseMatrix)).toBeCloseTo(0.0);

    expect(SparseMatrix.near(SparseMatrix.identity(10).inverse(), SparseMatrix.identity(10))).toBeTruthy();
    let permutationMatrix = new PermutationMatrix([1, 6, 8, 2, 5, 4, 9, 3, 0, 7], true);
    assert(permutationMatrix.isValid(), "Invalid permutation");
    expect(SparseMatrix.near(permutationMatrix.toSparseMatrix().inverse(), permutationMatrix.inverse().toSparseMatrix()));
    //expect(Math.abs(permutationMatrix.toSparseMatrix().determinant())).toBeCloseTo(1.0);
    //expect(nonSingularSparseMatrix.determinant()).not.toBeCloseTo(0.0);
    expect(SparseMatrix.near(SparseMatrix.mul(nonSingularSparseMatrix.inverse(), nonSingularSparseMatrix), SparseMatrix.identity(nonSingularDenseMatrix.width())));
    // test trivial singular matrix
    //expect(singularSparseTrivialMatrix.determinant()).toBeCloseTo(0.0);
    //expect(singularSparseNonTrivialMatrix.determinant()).toBeCloseTo(0.0);
});