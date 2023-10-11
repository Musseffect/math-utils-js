import Matrix from "../denseMatrix";
import mat2 from "../mat2";
import mat3 from "../mat3";
import mat4 from "../mat4";
import Triplet from "../triplet";
import { Tolerance, SmallTolerance, assert } from "../utils";
import vec2 from "../vec2";
import vec3 from "../vec3";
import vec4 from "../vec4";
import Vector from "../vector";


const singularTrivialMatrixTriplets: Triplet[] = [{ row: 1, column: 1, value: 1 },
{ row: 2, column: 2, value: 1 }, { row: 0, column: 0, value: 2 }, { row: 3, column: 3, value: 1 }, { row: 5, column: 5, value: 1 }];
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
const singularDenseNonTrivialMatrix = Matrix.fromTriplets(singularNonTrivialMatrixTriplets, 6, 6);

const nonSingularMatrixTriplets: Triplet[] = [{ row: 0, column: 0, value: 1 }, { row: 0, column: 2, value: 3 }, { row: 0, column: 4, value: 1 },
{ row: 1, column: 0, value: 2 }, { row: 1, column: 1, value: 2 }, { row: 1, column: 5, value: 2 }, { row: 2, column: 1, value: 2 }, { row: 2, column: 3, value: 1 }
    , { row: 2, column: 5, value: 3 }, { row: 3, column: 4, value: 4 }, { row: 4, column: 0, value: 3 }, { row: 4, column: 1, value: -3 }, { row: 4, column: 2, value: -1 }
    , { row: 4, column: 4, value: 2 }, { row: 5, column: 1, value: 2 }, { row: 5, column: 3, value: 2 }, { row: 5, column: 5, value: 4 }];
const nonSingularDenseMatrix = Matrix.fromTriplets(nonSingularMatrixTriplets, 6, 6);
const nonSingularMatrixDeterminant = 144.0;

test('Matrix operations', () => {
    expect(mat2.identity().determinant()).toBeCloseTo(1);
    expect(mat3.identity().determinant()).toBeCloseTo(1);
    expect(mat4.identity().determinant()).toBeCloseTo(1);
    expect(Matrix.identity(4).determinantNaive()).toBeCloseTo(1);

    expect(mat2.empty().determinant()).toBeCloseTo(0);
    expect(mat3.empty().determinant()).toBeCloseTo(0);
    expect(mat4.empty().determinant()).toBeCloseTo(0);
    expect(Matrix.empty(4, 4).determinantNaive()).toBeCloseTo(0);

    expect(Matrix.identity(2).determinantNaive()).toBeCloseTo(1);
    expect(Matrix.identity(3).determinantNaive()).toBeCloseTo(1);
    expect((new Matrix([1, 2, 3, 4], 2, 2)).determinantNaive()).toBeCloseTo(-2);

    const matT = new Matrix([3, 2, 4, 1, 2, 3, 1, 5, 2], 3, 3);
    expect(matT.determinantNaive()).toBeCloseTo(-19);

    const mat = new Matrix([4, 2, 0, 2, 0,
        2, 2, 0, 0, 0,
        3, 2, 3, 1, 0,
        0, 0, 0, 0, 4,
        0, -3, -1, 0, 2], 5, 5);
    console.log(mat.toString());
    expect(mat.determinantNaive()).toBeCloseTo(144);
    const matN = new Matrix([0, 2, -1, 0,
        2, 0, 0, 0, 3, 0, 3, 1, 0, 4, 0, 0], 4, 4);
    expect(matN.determinantNaive()).toBeCloseTo(8);
    expect(singularDenseTrivialMatrix.determinantNaive()).toBeCloseTo(0);
    expect(singularDenseNonTrivialMatrix.determinantNaive()).toBeCloseTo(0);
    expect(nonSingularDenseMatrix.determinantNaive()).toBeCloseTo(nonSingularMatrixDeterminant);

    let m4 = new mat4(
        1, 3, 2, 4,
        6, 8, 3, -2,
        -5, 3, 2, 1,
        3, 4, 5, 2);
    expect(Matrix.near(m4.transpose().transpose(), m4, Tolerance)).toBeTruthy();

    expect(Matrix.near(mat4.identity(), mat4.mul(m4, m4.inverse()), Tolerance)).toBeTruthy();

    let m3a = new mat3(1, 3, 2, 4, 6, 8, 3, -2, -5);
    let m3b = new mat3(4, -2, -1, -31, 21, -4, 51, -13, 10);
    let m3c = new mat3(13, 35, 7, 238, 14, 52, -181, 17, -45);
    expect(Matrix.near(m3a, mat3.mul(m3a, mat3.identity()), Tolerance)).toBeTruthy();
    expect(Matrix.near(m3c, mat3.mul(m3a, m3b), Tolerance)).toBeTruthy();
    expect(Matrix.near(mat3.identity(), mat3.mul(m3a, m3a.inverse()), Tolerance)).toBeTruthy();

    //let m2 = new mat2(-13.1, 0.6, -1.7, 2.3);
    //expect(mat2.near(mat2.identity(), mat2.mul(m2, m2.inverse()), Tolerance)).toBeTruthy();

});

test('Vector operations', () => {
    let a3D = new vec3(1., 3., 2.);
    let b3D = new vec3(2., -1., -4.);
    expect(a3D.l1norm()).toBeCloseTo(6, Tolerance);
    expect(a3D.l2norm()).toBeCloseTo(Math.sqrt(vec3.dot(a3D, a3D)), Tolerance);
    expect(a3D.lInfnorm()).toBeCloseTo(3, Tolerance);
    expect(vec3.dot(a3D, vec3.cross(a3D, b3D))).toBeCloseTo(0.0);
    expect(vec3.near(vec3.add(a3D, b3D), new vec3(3, 2, -2), Tolerance)).toBeTruthy();
    expect(vec3.near(vec3.mul(a3D, b3D), new vec3(2, -3, -8), Tolerance)).toBeTruthy();
    expect(vec3.near(vec3.sub(a3D, b3D), new vec3(-1, 4, 6), Tolerance)).toBeTruthy();
    expect(vec3.near(vec3.div(a3D, b3D), new vec3(0.5, -3, -0.5), Tolerance)).toBeTruthy();
});

test('Matrix-vector operations', () => {
    let p2D = new vec2(2., -3.);
    let p3D = new vec3(1.6, 0.5, -0.2);
    let p4D = new vec4(-1., 0.3, 2., -4.);
    let m2D = new mat2(
        2., 0.3,
        7., -13.
    );
    let m3D = new mat3(
        0.3, 1.3, -60.3,
        -0.9, 0.51, 3.1,
        10.23, 1.2, 3.4
    );
    let m4D = new mat4(
        0.3, 1.3, -60.3, 0.12,
        -0.9, 0.51, 3.1, 14,
        10.23, 1.2, 3.4, 1.7,
        0, -6.75, 1, -0.8
    );
    let m2DInv = m2D.inverse();
    let m3DInv = m3D.inverse();
    let m4DInv = m4D.inverse();
    expect(vec2.near(p2D, m2DInv.postMulVec(m2D.postMulVec(p2D))));
    expect(vec2.near(p2D, m2DInv.preMulVec(m2D.preMulVec(p2D))));
    expect(vec3.near(p3D, m3DInv.postMulVec(m3D.postMulVec(p3D))));
    expect(vec3.near(p3D, m3DInv.preMulVec(m3D.preMulVec(p3D))));
    expect(vec4.near(p4D, m4DInv.postMulVec(m4D.postMulVec(p4D))));
    expect(vec4.near(p4D, m4DInv.preMulVec(m4D.preMulVec(p4D))));

    /*m4D = new mat4(
        0, 1, 0, 1,
        1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 1, 0
    );*/
    let mat = m4D.toMatrix();
    let rhs = p4D.toVector();
    expect(Vector.near(Matrix.postMulVec(mat, rhs), m4D.postMulVec(p4D).toVector(), SmallTolerance)).toBeTruthy();
    expect(Vector.near(Matrix.preMulVec(mat, rhs), m4D.preMulVec(p4D).toVector(), SmallTolerance)).toBeTruthy();

    //expect(new Matrix([1, 2, 3, 4], 2, 2).determinantNaive()).toBeCloseTo(-2.0);

    expect(new Matrix([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3).determinantNaive()).toBeCloseTo(0.0);
    //expect(mat.determinantNaive()).toBeCloseTo(m4D.determinant());
    //console.log(`InverseNaive ${mat.inverseNaive().toString()}`);
    console.log(`Inverse ${m4D.inverse().toString()}`);
    //expect(Matrix.near(mat.inverseNaive(), m4D.inverse().toMatrix(), SmallTolerance)).toBeTruthy();

    //expect(Matrix.near(Matrix.mul(mat, mat.inverseNaive()), Matrix.identity(4), SmallTolerance)).toBeTruthy();
    expect(Vector.near(m4D.inverse().postMulVec(p4D).toVector(), Matrix.solve(mat.clone(), rhs), SmallTolerance)).toBeTruthy();
    expect(Matrix.near(mat.transpose(), m4D.transpose().toMatrix(), SmallTolerance)).toBeTruthy();
});
test("General dense matrix", () => {
    // test nonsingular matrices
    expect(Matrix.near(Matrix.identity(4).inverseNaive(), Matrix.identity(4))).toBeTruthy();
    //expect(nonSingularDenseMatrix.determinant()).not.toBeCloseTo(0.0);
    //expect(Matrix.near(Matrix.mul(nonSingularDenseMatrix.inverse(), nonSingularDenseMatrix), Matrix.identity(nonSingularDenseMatrix.width())));
    // test trivial singular matrix
    //expect(singularDenseTrivialMatrix.determinant()).toBeCloseTo(0.0);
    //expect(singularDenseNonTrivialMatrix.determinant()).toBeCloseTo(0.0);

    // test non-trivial singular matrix
});

test.skip('Triangular matrix', () => {
    // todo
    throw new Error("Not implemented");
});

