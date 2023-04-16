import axisAngle from "./axisAngle";
import mat3 from "./mat3";
import mat4 from "./mat4";
import quat from "./quat";
import { assert, Epsilon, near, radians, SmallEpsilon } from "./utils";
import transform3D from "./transform3D";
import transform2D from "./transform2D";
import complex from "./complex";
import vec2 from "./vec2";
import vec3 from "./vec3";
import mat2 from "./mat2";
import vec4 from "./vec4";
import Matrix from "./denseMatrix";
import Vector from "./vector";
import sparseVector from "./sparseVector";
import PermutationMatrix from "./permutationMatrix";
import SparseMatrix from "./sparseMatrix";
import Triplet from "./triplet";
import { backwardDifference, centralDifference, forwardDifference, secondOrderDifference } from "./numericalDifferentiation";

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


test.only('Quaternion basic operations', () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let q1 = quat.fromAxisAngle(axisAngleRotation);
    expect(quat.near(q1.mul(quat.inverse(q1)), quat.identity())).toBeTruthy();
    let q2 = quat.fromComponents(1.0, 2.0, 0.3, 0.157);
    let qDiv = q1.div(q2);
    let qAdd = q1.add(q2);
    let qSub = q1.sub(q2);
    let qMul = q1.mul(q2);
    expect(quat.near(q1, qDiv.mul(q2))).toBeTruthy();
    expect(quat.near(q1, qAdd.sub(q2))).toBeTruthy();
    expect(quat.near(q1, qSub.add(q2))).toBeTruthy();
    expect(quat.near(q1, qMul.div(q2))).toBeTruthy();
});

test('Complex operations', () => {
    let z1 = new complex(1, 2);
    let z2 = complex.polar(2.0, 1.1);
    expect(z2.length()).toBeCloseTo(2.0);
    expect(z2.arg()).toBeCloseTo(1.1);
    let zDiv = complex.cDiv(z1, z2);
    let zMul = complex.cMul(z1, z2);
    let zExp = complex.cExp(z1);
    let zLog = complex.cLog(z1);
    let zPow = complex.cPow(z1, z2);
    expect(vec2.near(z1, complex.cMul(zDiv, z2))).toBeTruthy();
    expect(vec2.near(z1, complex.cDiv(zMul, z2))).toBeTruthy();
    expect(vec2.near(z1, complex.cLog(zExp))).toBeTruthy();
    expect(vec2.near(z1, complex.cExp(zLog))).toBeTruthy();
    expect(vec2.near(z2, complex.cDiv(complex.cLog(zPow), complex.cLog(z1)))).toBeTruthy();
});

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

    expect(Matrix.lInfDistance(singularSparseTrivialMatrix.toDense(), singularDenseTrivialMatrix)).toBeCloseTo(0.0);
    expect(Matrix.lInfDistance(singularSparseNonTrivialMatrix.toDense(), singularDenseNonTrivialMatrix)).toBeCloseTo(0.0);
    expect(Matrix.lInfDistance(nonSingularSparseMatrix.toDense(), nonSingularDenseMatrix)).toBeCloseTo(0.0);

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
    expect(Matrix.near(m4.transpose().transpose(), m4, Epsilon)).toBeTruthy();

    expect(Matrix.near(mat4.identity(), mat4.mul(m4, m4.inverse()), Epsilon)).toBeTruthy();

    let m3a = new mat3(1, 3, 2, 4, 6, 8, 3, -2, -5);
    let m3b = new mat3(4, -2, -1, -31, 21, -4, 51, -13, 10);
    let m3c = new mat3(13, 35, 7, 238, 14, 52, -181, 17, -45);
    expect(Matrix.near(m3a, mat3.mul(m3a, mat3.identity()), Epsilon)).toBeTruthy();
    expect(Matrix.near(m3c, mat3.mul(m3a, m3b), Epsilon)).toBeTruthy();
    expect(Matrix.near(mat3.identity(), mat3.mul(m3a, m3a.inverse()), Epsilon)).toBeTruthy();

    //let m2 = new mat2(-13.1, 0.6, -1.7, 2.3);
    //expect(mat2.near(mat2.identity(), mat2.mul(m2, m2.inverse()), Epsilon)).toBeTruthy();

});

test('Vector operations', () => {
    let a3D = new vec3(1., 3., 2.);
    let b3D = new vec3(2., -1., -4.);
    expect(a3D.l1norm()).toBeCloseTo(6, Epsilon);
    expect(a3D.l2norm()).toBeCloseTo(Math.sqrt(vec3.dot(a3D, a3D)), Epsilon);
    expect(a3D.lInfnorm()).toBeCloseTo(3, Epsilon);
    expect(vec3.dot(a3D, vec3.cross(a3D, b3D))).toBeCloseTo(0.0);
    expect(vec3.near(vec3.add(a3D, b3D), new vec3(3, 2, -2), Epsilon)).toBeTruthy();
    expect(vec3.near(vec3.mul(a3D, b3D), new vec3(2, -3, -8), Epsilon)).toBeTruthy();
    expect(vec3.near(vec3.sub(a3D, b3D), new vec3(-1, 4, 6), Epsilon)).toBeTruthy();
    expect(vec3.near(vec3.div(a3D, b3D), new vec3(0.5, -3, -0.5), Epsilon)).toBeTruthy();
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
    expect(Vector.near(Matrix.postMulVec(mat, rhs), m4D.postMulVec(p4D).toVector(), SmallEpsilon)).toBeTruthy();
    expect(Vector.near(Matrix.preMulVec(mat, rhs), m4D.preMulVec(p4D).toVector(), SmallEpsilon)).toBeTruthy();

    //expect(new Matrix([1, 2, 3, 4], 2, 2).determinantNaive()).toBeCloseTo(-2.0);

    expect(new Matrix([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3).determinantNaive()).toBeCloseTo(0.0);
    //expect(mat.determinantNaive()).toBeCloseTo(m4D.determinant());
    //console.log(`InverseNaive ${mat.inverseNaive().toString()}`);
    console.log(`Inverse ${m4D.inverse().toString()}`);
    //expect(Matrix.near(mat.inverseNaive(), m4D.inverse().toMatrix(), SmallEpsilon)).toBeTruthy();

    //expect(Matrix.near(Matrix.mul(mat, mat.inverseNaive()), Matrix.identity(4), SmallEpsilon)).toBeTruthy();
    expect(Vector.near(m4D.inverse().postMulVec(p4D).toVector(), Matrix.solve(mat.clone(), rhs), SmallEpsilon)).toBeTruthy();
    expect(Matrix.near(mat.transpose(), m4D.transpose().toMatrix(), SmallEpsilon)).toBeTruthy();
});

test("Rotation conversions", () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));

    let quatRotation = quat.fromAxisAngle(axisAngleRotation);
    let matRotation = mat3.fromAxisAngle(axisAngleRotation);
    expect(matRotation.determinant()).toBeCloseTo(1);
    expect(Matrix.near(quatRotation.toMat3(), matRotation, SmallEpsilon)).toBeTruthy();
    expect(quat.near(quatRotation, matRotation.toQuat(), SmallEpsilon)).toBeTruthy();
    expect(axisAngle.near(quatRotation.toAxisAngle(), axisAngleRotation, SmallEpsilon)).toBeTruthy();
    expect(axisAngle.near(matRotation.toAxisAngle(), axisAngleRotation, SmallEpsilon)).toBeTruthy();

    let point = new vec3(0.3, -0.5, -0.2);
    let eulerRotation = { yaw: radians(10), pitch: radians(25), roll: radians(62) };
    matRotation = mat3.fromEulerAngles(eulerRotation.yaw, eulerRotation.pitch, eulerRotation.roll);
    quatRotation = quat.fromEulerAngles(eulerRotation.yaw, eulerRotation.pitch, eulerRotation.roll);
    let matRotationEuler = mat3.yaw(eulerRotation.yaw).mulSelf(mat3.pitch(eulerRotation.pitch)).mulSelf(mat3.roll(eulerRotation.roll));
    expect(vec3.near(quatRotation.toMat3().transformPoint3D(point), matRotation.transformPoint3D(point), SmallEpsilon)).toBeTruthy();
    expect(vec3.near(quatRotation.rotate(point), matRotation.toQuat().rotate(point), SmallEpsilon)).toBeTruthy();
    expect(vec3.near(matRotation.transformPoint3D(point), matRotationEuler.transformPoint3D(point), SmallEpsilon)).toBeTruthy();
    expect(Matrix.near(matRotation, matRotationEuler, SmallEpsilon)).toBeTruthy();

    matRotation = mat3.fromEulerAngles(eulerRotation.yaw, eulerRotation.pitch, eulerRotation.roll);
    let ea = matRotation.toEulerAngles();
    let eaMat = mat3.fromEulerAngles(ea.x, ea.y, ea.z);

    expect(vec3.near(eaMat.transformPoint3D(point), matRotation.transformPoint3D(point), SmallEpsilon)).toBeTruthy();

    matRotation = mat3.fromEulerAngles(radians(180), 0, 0);
    expect(Math.abs(vec3.dot(matRotation.axis(), new vec3(0, 1, 0)))).toBeCloseTo(1, 4);
});

test("Rotations", () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let trans = new transform3D(new vec3(-2, 3, 4.), quat.fromAxisAngle(axisAngleRotation), new vec3(1.2, -0.2, 0.4));
    let point = new vec3(0.3, -0.5, -0.2);
    let vector = new vec3(0.3, -0.5, -0.2);
    let affineTransform = trans.toAffine();
    expect(vec3.near(affineTransform.transformPoint3D(point), trans.transformPoint3D(point), Epsilon)).toBeTruthy();
    expect(vec3.near(affineTransform.transformVector3D(vector), trans.transformVector3D(vector), Epsilon)).toBeTruthy();

    let quatRotation = quat.fromAxisAngle(axisAngleRotation);
    let matRotation = mat3.fromAxisAngle(axisAngleRotation);
    expect(vec3.near(quatRotation.rotate(point), matRotation.transformPoint3D(point), Epsilon)).toBeTruthy();
    expect(vec3.near(quatRotation.rotate(point), axisAngleRotation.rotate(point), Epsilon)).toBeTruthy();

    let eulerRotation = new vec3(radians(30), radians(55), radians(72));
    quatRotation = quat.fromEulerAngles(eulerRotation.y, eulerRotation.x, eulerRotation.z);
    matRotation = mat3.fromEulerAngles(eulerRotation.y, eulerRotation.x, eulerRotation.z);
    axisAngleRotation = quatRotation.toAxisAngle();
    expect(vec3.near(quatRotation.rotate(point), matRotation.transformPoint3D(point), Epsilon)).toBeTruthy();
    expect(vec3.near(quatRotation.rotate(point), eulerRotation.rotateEuler(point), Epsilon)).toBeTruthy();
    expect(vec3.near(axisAngleRotation.rotate(point), eulerRotation.rotateEuler(point), Epsilon)).toBeTruthy();

    axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    matRotation = mat3.fromAxisAngle(axisAngleRotation);
    expect(Matrix.near(matRotation.inverse(), matRotation.transpose(), Epsilon)).toBeTruthy();
});

test("Transformations", () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let translation3D = new vec3(-2, 3, 4.);
    let rotation3D = quat.fromAxisAngle(axisAngleRotation);
    let scale3D = new vec3(1.2, -0.2, 0.4);
    let t3D = mat4.fromTranslation3D(translation3D);
    let r3D = mat4.fromRotation3D(rotation3D);
    let s3D = mat4.fromScale3D(scale3D);
    let affine3D = mat4.fromTRS(translation3D, rotation3D, scale3D);
    let inverseAffine3D = affine3D.inverse();
    let transformation3D = new transform3D(translation3D, rotation3D, scale3D);
    let point3D = new vec3(0.6, -1.2, 1.6);
    expect(Matrix.near(affine3D, t3D.mul(r3D).mul(s3D), Epsilon)).toBeTruthy();
    expect(vec3.near(transformation3D.transformPoint3D(point3D), affine3D.transformPoint3D(point3D), Epsilon)).toBeTruthy();
    expect(vec3.near(transformation3D.transformVector3D(point3D), affine3D.transformVector3D(point3D), Epsilon)).toBeTruthy();
    expect(vec3.near(transformation3D.transformInversePoint3D(point3D), inverseAffine3D.transformPoint3D(point3D), Epsilon)).toBeTruthy();
    expect(vec3.near(transformation3D.transformInverseVector3D(point3D), inverseAffine3D.transformVector3D(point3D), Epsilon)).toBeTruthy();
    expect(transformation3D.hasUniformScale()).toBeFalsy();


    let rigidTransform3D = new transform3D(translation3D, rotation3D, new vec3(1., 1., 1.));
    let inverseRigidTransform3D = rigidTransform3D.inverse();
    expect(vec3.near(point3D, affine3D.transformPoint3D(inverseAffine3D.transformPoint3D(point3D)), Epsilon)).toBeTruthy();
    expect(vec3.near(rigidTransform3D.transformInversePoint3D(point3D),
        inverseRigidTransform3D.transformPoint3D(point3D), SmallEpsilon)).toBeTruthy();

    let translation2D = new vec2(1.0, -2.0);
    let rotation2D = radians(30);
    let scale2D = new vec2(-1.1, 0.1);
    let t2D = mat3.fromTranslation2D(translation2D);
    let r2D = mat3.fromRotation2D(rotation2D);
    let s2D = mat3.fromScale2D(scale2D);
    let affine2D = mat3.fromTRS(translation2D, rotation2D, scale2D);
    let inverseAffine2D = affine2D.inverse();
    let transformation2D = new transform2D(translation2D, rotation2D, scale2D);
    let point2D = new vec2(3.6, 1.6);
    expect(Matrix.near(affine2D, t2D.mul(r2D).mul(s2D), Epsilon)).toBeTruthy();
    expect(vec2.near(transformation2D.transformPoint2D(point2D), affine2D.transformPoint2D(point2D), Epsilon)).toBeTruthy();
    expect(vec2.near(transformation2D.transformVector2D(point2D), affine2D.transformVector2D(point2D), Epsilon)).toBeTruthy();

    expect(vec2.near(point2D, affine2D.transformPoint2D(inverseAffine2D.transformPoint2D(point2D)), Epsilon)).toBeTruthy();

    expect(vec2.near(transformation2D.transformInversePoint2D(point2D), inverseAffine2D.transformPoint2D(point2D), Epsilon)).toBeTruthy();
    expect(vec2.near(transformation2D.transformInverseVector2D(point2D), inverseAffine2D.transformVector2D(point2D), Epsilon)).toBeTruthy();
    expect(transformation2D.hasUniformScale()).toBeFalsy();

    let rigidTransform2D = new transform2D(translation2D, rotation2D, new vec2(-1., -1.));
    let inverseRigidTransform2D = rigidTransform2D.inverse();
    expect(vec2.near(rigidTransform2D.transformInversePoint2D(point2D),
        inverseRigidTransform2D.transformPoint2D(point2D), Epsilon)).toBeTruthy();
});


test("Transform components", () => {
    let point3D = new vec3(0.3, -0.5, -0.2);
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let translation3D = new vec3(-2, 3., 4.);
    let rotation3D = quat.fromAxisAngle(axisAngleRotation);
    let scale3D = new vec3(1.2, -0.2, 0.4);
    let t = mat4.fromTranslation3D(translation3D);
    let r = mat4.fromRotation3D(rotation3D);
    let s = mat4.fromScale3D(scale3D);
    let trs = mat4.fromTRS(translation3D, rotation3D, scale3D);
    let pS1 = vec3.mul(point3D, scale3D);
    let pS2 = s.transformPoint3D(point3D);
    let pRS1 = rotation3D.rotate(pS1);
    let pRS2 = r.transformPoint3D(pS2);
    let pTRS1 = vec3.add(translation3D, pRS1);
    let pTRS2 = t.transformPoint3D(pRS2);
    expect(vec3.near(pS1, pS2, Epsilon)).toBeTruthy();
    expect(vec3.near(pRS1, pRS2, Epsilon)).toBeTruthy();
    expect(vec3.near(pTRS1, pTRS2, Epsilon)).toBeTruthy();
    expect(vec3.near(pTRS1, trs.transformPoint3D(point3D), Epsilon)).toBeTruthy();

    let point2D = new vec2(3.0, 13.0);
    let translation2D = new vec2(-3.0, 0.1);
    let rotation2D = radians(30.0);
    let scale2D = new vec2(0.1, -2.3);
    let t2D = mat3.fromTranslation2D(translation2D);
    let r2D = mat3.fromRotation2D(rotation2D);
    let s2D = mat3.fromScale2D(scale2D);
    let trs_2D = mat3.fromTRS(translation2D, rotation2D, scale2D);
    let pS1_2D = vec2.mul(point2D, scale2D);
    let pS2_2D = s2D.transformPoint2D(point2D);
    let pRS1_2D = vec2.rotate2D(pS1_2D, rotation2D);
    let pRS2_2D = r2D.transformPoint2D(pS2_2D);
    let pTRS1_2D = vec2.add(translation2D, pRS1_2D);
    let pTRS2_2D = t2D.transformPoint2D(pRS2_2D);
    expect(vec2.near(pS1_2D, pS1_2D, SmallEpsilon)).toBeTruthy();
    expect(vec2.near(pRS1_2D, pRS2_2D, SmallEpsilon)).toBeTruthy();
    expect(vec2.near(pTRS1_2D, pTRS2_2D, SmallEpsilon)).toBeTruthy();
    expect(vec2.near(pTRS1_2D, trs_2D.transformPoint2D(point2D), SmallEpsilon)).toBeTruthy();
});



test("General dense matrix", () => {

    // test nonsingular matrices
    expect(Matrix.near(Matrix.identity(4).inverseNaive(), Matrix.identity(4))).toBeTruthy();
    let permutationMatrix = new PermutationMatrix([1, 6, 8, 2, 5, 4, 9, 3, 0, 7], true);
    assert(permutationMatrix.isValid(), "Invalid permutation");
    //expect(Matrix.near(permutationMatrix.toMatrix().inverse(), permutationMatrix.inverse().toMatrix()));
    //expect(Math.abs(permutationMatrix.toMatrix().determinant())).toBeCloseTo(1.0);
    //expect(nonSingularDenseMatrix.determinant()).not.toBeCloseTo(0.0);
    //expect(Matrix.near(Matrix.mul(nonSingularDenseMatrix.inverse(), nonSingularDenseMatrix), Matrix.identity(nonSingularDenseMatrix.width())));
    // test trivial singular matrix
    //expect(singularDenseTrivialMatrix.determinant()).toBeCloseTo(0.0);
    //expect(singularDenseNonTrivialMatrix.determinant()).toBeCloseTo(0.0);

    // test non-trivial singular matrix
});

test.skip('Triangular matrix', () => { });


test("Sparse vector", () => {
    let dense: Vector = new Vector([0, -1, 2, -1e-8, 0.0, 3, -5, .0]);
    let v: sparseVector = sparseVector.fromVector(dense.data, SmallEpsilon);
    expect(Vector.near(v.toDense(), dense, SmallEpsilon)).toBeTruthy();
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

    v = sparseVector.fromVector(dense.data, SmallEpsilon);
    let dense2: Vector = new Vector([3, 0.1, 0.0, 1e-7, 0.3, -2.0, -6.0, 7.0]);
    let v2 = sparseVector.fromVector(dense2.data, SmallEpsilon);
    expect(Vector.near(v2.toDense(), dense2, SmallEpsilon));
    console.log(sparseVector.add(v, v2).toDense().toString());
    console.log(Vector.add(dense, dense2).toString());
    expect(Vector.near(sparseVector.add(v, v2).toDense(), Vector.add(dense, dense2))).toBeTruthy();
    expect(Vector.near(sparseVector.sub(v, v2).toDense(), Vector.sub(dense, dense2))).toBeTruthy();
    expect(Vector.near(sparseVector.mul(v, v2).toDense(), Vector.mul(dense, dense2))).toBeTruthy();
    expect(sparseVector.dot(v, v2),).toBeCloseTo(Vector.dot(dense, dense2));

    //expect(Vector.near(sparseVector.div(v, v2).toDense(), Vector.div(dense, dense2))).toBeTruthy();
});

test.skip("General sparse matrix", () => {

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


test("Differentiation", () => {
    let testFunc =
    {
        f: (x: Vector) => {
            return x.get(0) * x.get(0) + Math.sin(3.0 * x.get(0) * x.get(1));
        },
        dfdx: (x: Vector): Vector => {
            const cos = Math.cos(3.0 * x.get(0) * x.get(1));
            return new Vector([2.0 * x.get(0) - 3.0 * x.get(1) * cos, -3.0 * x.get(0) * cos]);
        },
        dfdxdy: (x: Vector): Matrix => {
            const sin = Math.sin(3.0 * x.get(0) * x.get(1));
            const cos = Math.cos(3.0 * x.get(0) * x.get(1));
            let result = Matrix.empty(2, 2);
            result.set(0, 0, 2.0 - 9.0 * x.get(1) * x.get(1) * sin);
            result.set(0, 1, -3.0 * cos - 9.0 * x.get(0) * x.get(1) * sin);
            result.set(1, 0, -3.0 * cos - 9.0 * x.get(0) * x.get(1) * sin);
            result.set(1, 1, -9.0 * x.get(0) * x.get(0) * sin);
            return result;
        }
    };
    let points: Vector[] = [
        new Vector([0, 0]), new Vector([-25, 12]), new Vector([-3, -3]), new Vector([0, 5]),
        new Vector([5, 0]), new Vector([14, 14])];
    const step = SmallEpsilon;
    const tolerance = Epsilon
    for (let point of points) {
        const dfdx = testFunc.dfdx(point);
        const dfdxdy = testFunc.dfdxdy(point);
        expect(Vector.sub(dfdx, forwardDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
        expect(Vector.sub(dfdx, backwardDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
        expect(Vector.sub(dfdx, centralDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
        expect(Matrix.sub(dfdxdy, secondOrderDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
    }
});