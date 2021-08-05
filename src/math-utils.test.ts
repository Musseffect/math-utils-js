import axisAngle from "./axisAngle";
import complex from "./complex";
import mat3 from "./mat3";
import mat4 from "./mat4";
import quat from "./quat";
import { Epsilon, near, radians, SmallEpsilon } from "./utils";
import transform from "./transform";
import vec2 from "./vec2";
import vec3 from "./vec3";
import mat2 from "./mat2";

test('Quaternion basic operations', () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let q1 = quat.fromAxisAngle(axisAngleRotation);
    expect(quat.near(q1.mul(q1.inverse()), quat.identity())).toBeTruthy();
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
    expect(mat4.identity().determinant()).toBeCloseTo(1);
    expect(mat3.identity().determinant()).toBeCloseTo(1);
    expect(mat2.identity().determinant()).toBeCloseTo(1);

    let m4 = new mat4(
        1, 3, 2, 4,
        6, 8, 3, -2,
        -5, 3, 2, 1,
        3, 4, 5, 2);
    expect(mat4.near(m4.transpose().transpose(), m4, Epsilon)).toBeTruthy();

    expect(mat4.near(mat4.identity(), mat4.mul(m4, m4.inverse()), Epsilon)).toBeTruthy();

    let m3a = new mat3(1, 3, 2, 4, 6, 8, 3, -2, -5);
    let m3b = new mat3(4, -2, -1, -31, 21, -4, 51, -13, 10);
    let m3c = new mat3(13, 35, 7, 238, 14, 52, -181, 17, -45);
    expect(mat3.near(m3a, mat3.mul(m3a, mat3.identity()), Epsilon)).toBeTruthy();
    expect(mat3.near(m3c, mat3.mul(m3a, m3b), Epsilon)).toBeTruthy();
});

test('Vector operations', () => {
    let a = new vec3(1., 3., 2.);
    let b = new vec3(2., -1., -4.);
    expect(a.l1norm()).toBeCloseTo(6, Epsilon);
    expect(a.l2norm()).toBeCloseTo(Math.sqrt(vec3.dot(a, a)), Epsilon);
    expect(a.lInfnorm()).toBeCloseTo(3, Epsilon);
    expect(vec3.dot(a, vec3.cross(a, b))).toBeCloseTo(0.0);
    expect(vec3.near(vec3.add(a, b), new vec3(3, 2, -2), Epsilon)).toBeTruthy();
    expect(vec3.near(vec3.mul(a, b), new vec3(2, -3, -8), Epsilon)).toBeTruthy();
    expect(vec3.near(vec3.sub(a, b), new vec3(-1, 4, 6), Epsilon)).toBeTruthy();
    expect(vec3.near(vec3.div(a, b), new vec3(0.5, -3, -0.5), Epsilon)).toBeTruthy();
});

test('Matrix-vector operations', () => {
});

test("Rotation conversions test", () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));

    let quatRotation = quat.fromAxisAngle(axisAngleRotation);
    let matRotation = mat3.fromAxisAngle(axisAngleRotation);
    let point = new vec3(0.3, -0.5, -0.2);
    expect(mat3.near(quatRotation.toMat3(), matRotation, SmallEpsilon)).toBeTruthy();
    expect(quat.near(quatRotation, matRotation.toQuat(), SmallEpsilon)).toBeTruthy();
    expect(axisAngle.near(quatRotation.toAxisAngle(), axisAngleRotation, SmallEpsilon)).toBeTruthy();
    expect(axisAngle.near(matRotation.toAxisAngle(), axisAngleRotation, SmallEpsilon)).toBeTruthy();

    let eulerRotation = new vec3(radians(30), radians(55), radians(72));
    quatRotation = quat.fromEulerAngles(eulerRotation.x, eulerRotation.y, eulerRotation.z);
    matRotation = mat3.fromEulerAngles(eulerRotation.x, eulerRotation.y, eulerRotation.z);
    expect(vec3.near(quatRotation.toMat3().transform(point), matRotation.transform(point), SmallEpsilon)).toBeTruthy();
    expect(vec3.near(quatRotation.rotate(point), matRotation.toQuat().rotate(point), SmallEpsilon)).toBeTruthy();
});

test("Rotations",() => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let trans = new transform(new vec3(-2, 3, 4.), quat.fromAxisAngle(axisAngleRotation), new vec3(1.2, -0.2, 0.4));
    let point = new vec3(0.3, -0.5, -0.2);
    let vector = new vec3(0.3, -0.5, -0.2);
    let affineTransform = trans.toAffine();
    expect(vec3.near(affineTransform.transformPoint(point), trans.transformPoint(point), Epsilon)).toBeTruthy();
    expect(vec3.near(affineTransform.transformVector(vector), trans.transformVector(vector), Epsilon)).toBeTruthy();

    let quatRotation = quat.fromAxisAngle(axisAngleRotation);
    let matRotation = mat3.fromAxisAngle(axisAngleRotation);
    expect(vec3.near(quatRotation.rotate(point), matRotation.transform(point), Epsilon)).toBeTruthy();
    expect(vec3.near(quatRotation.rotate(point), axisAngleRotation.rotate(point), Epsilon)).toBeTruthy();   

    let eulerRotation = new vec3(radians(30), radians(55), radians(72));
    quatRotation = quat.fromEulerAngles(eulerRotation.x, eulerRotation.y, eulerRotation.z);
    matRotation = mat3.fromEulerAngles(eulerRotation.x, eulerRotation.y, eulerRotation.z);
    axisAngleRotation = quatRotation.toAxisAngle();
    expect(vec3.near(quatRotation.rotate(point), matRotation.transform(point), Epsilon)).toBeTruthy();
    expect(vec3.near(quatRotation.rotate(point), eulerRotation.rotateEuler(point), Epsilon)).toBeTruthy();
    expect(vec3.near(axisAngleRotation.rotate(point), eulerRotation.rotateEuler(point), Epsilon)).toBeTruthy();

    axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    matRotation = mat3.fromAxisAngle(axisAngleRotation);
    expect(mat3.near(matRotation.inverse(), matRotation.transpose(), Epsilon)).toBeTruthy();
});

test("Transformations", () => {
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let translation = new vec3(-2, 3, 4.);
    let rotation = quat.fromAxisAngle(axisAngleRotation);
    let scale = new vec3(1.2, -0.2, 0.4);
    let t = mat4.translation(translation);
    let r = mat4.rotation(rotation);
    let s = mat4.fromScale(scale);
    let trs = mat4.fromTRS(translation, rotation, scale);
    expect(mat4.near(trs, t.mul(r).mul(s), Epsilon)).toBeTruthy();
});


test("Transform components", ()=>{
    let point = new vec3(0.3, -0.5, -0.2);
    let axisAngleRotation = new axisAngle(new vec3(1., 2., -3.), radians(70));
    let translation = new vec3(-2, 3, 4.);
    let rotation = quat.fromAxisAngle(axisAngleRotation);
    let scale = new vec3(1.2, -0.2, 0.4);
    let t = mat4.translation(translation);
    let r = mat4.rotation(rotation);
    let s = mat4.fromScale(scale);
    let trs = mat4.fromTRS(translation, rotation, scale);
    let pS1 = vec3.mul(point, scale);
    let pS2 = s.transformPoint(point);
    let pRS1 = rotation.rotate(pS1);
    let pRS2 = r.transformPoint(pS2);
    let pTRS1 = vec3.add(translation, pRS1);
    let pTRS2 = t.transformPoint(pRS2);
    expect(vec3.near(pS1, pS2, Epsilon)).toBeTruthy();
    expect(vec3.near(pRS1, pRS2, Epsilon)).toBeTruthy();
    expect(vec3.near(pTRS1, pTRS2, Epsilon)).toBeTruthy();
    expect(vec3.near(pTRS1, trs.transformPoint(point), Epsilon)).toBeTruthy();

});