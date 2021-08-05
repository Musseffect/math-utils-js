import axisAngle from "./axisAngle";
import mat2 from "./mat2";
import mat4 from "./mat4";
import quat from "./quat";
import { assert, determinant2x2, determinant3x3, Epsilon, near, SmallEpsilon } from "./utils";
import vec2 from "./vec2";
import vec3 from "./vec3";


export default class mat3 {
    data: number[][];
    constructor(m11: number, m12: number, m13: number,
        m21: number, m22: number, m23: number,
        m31: number, m32: number, m33: number) {
        this.data = [
            [m11, m12, m13],
            [m21, m22, m23],
            [m31, m32, m33]
        ];
    }
    clone(): mat3 {
        return new mat3(
            this.get(0, 0), this.get(0, 1), this.get(0, 2),
            this.get(1, 0), this.get(1, 1), this.get(1, 2),
            this.get(2, 0), this.get(2, 1), this.get(2, 2),
        );
    }
    static near(a: mat3, b: mat3, threshold?: number): boolean {
        if (!threshold) {
            threshold = Epsilon;
        }
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                if (Math.abs(a.get(i, j) - b.get(i, j)) > threshold)
                    return false;
            }
        }
        return true;
    }
    static empty(): mat3 {
        return new mat3(
            0, 0, 0,
            0, 0, 0,
            0, 0, 0);
    }
    static identity(): mat3 {
        return new mat3(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1);
    }
    static fromAxisAngle(aa: axisAngle): mat3 {
        let c = Math.cos(aa.angle);
        let s = Math.sin(aa.angle);
        let x = aa.axis.x;
        let y = aa.axis.y;
        let z = aa.axis.z;
        return new mat3(
            c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s,
            y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s,
            z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c));


        return new mat3(
            c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s,
            y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s,
            z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c));
    }
    trace(): number {
        return this.get(0, 0) + this.get(1, 1) + this.get(2, 2);
    }
    toAxisAngle(): axisAngle {
        let trace = this.trace();
        if (near(trace, 3, SmallEpsilon)) {
            let axis = new vec3(
                this.get(2, 1) - this.get(1, 2),
                this.get(0, 2) - this.get(2, 0),
                this.get(1, 0) - this.get(0, 1)
            );
            let angle = Math.acos((trace - 1) / 2);
            return new axisAngle(axis, angle);
        } else {
            let axis = new vec3(
                this.get(2, 1) - this.get(1, 2),
                this.get(0, 2) - this.get(2, 0),
                this.get(1, 0) - this.get(0, 1)
            );
            let length = axis.l2norm();
            let angle = Math.atan2(length, trace - 1);
            if (Math.abs(angle) < SmallEpsilon)
                axis = new vec3(1., 0., 0.);
            return new axisAngle(axis, angle);
        }
    }
    toEulerAngles(): vec3 {
        throw new Error("Not implemented");
    }
    //TODO: rewrite
    //y - up
    //x - left
    //z - dir 
    static yaw(yaw: number): mat3 {
        let cy = Math.cos(yaw);
        let sy = Math.cos(yaw);
        let m11 = 1;
        let m12 = 0;
        let m13 = 0;
        let m21 = 0;
        let m22 = cy;
        let m23 = -sy;
        let m31 = 0;
        let m32 = sy;
        let m33 = cy;
        return new mat3(
            m11, m12, m13,
            m21, m22, m23,
            m31, m32, m33);
    }
    static pitch(pitch: number): mat3 {
        let cp = Math.cos(pitch);
        let sp = Math.cos(pitch);
        let m11 = cp;
        let m12 = 0;
        let m13 = sp;
        let m21 = 0;
        let m22 = 1;
        let m23 = 0;
        let m31 = -sp;
        let m32 = 0;
        let m33 = cp;
        return new mat3(
            m11, m12, m13,
            m21, m22, m23,
            m31, m32, m33);
    }
    static roll(roll: number): mat3 {
        let cr = Math.cos(roll);
        let sr = Math.cos(roll);
        let m11 = cr;
        let m12 = -sr;
        let m13 = 0;
        let m21 = sr;
        let m22 = cr;
        let m23 = 0;
        let m31 = 0;
        let m32 = 0;
        let m33 = 1;
        return new mat3(
            m11, m12, m13,
            m21, m22, m23,
            m31, m32, m33);
    }
    static fromEulerAngles(yaw: number, pitch: number, roll: number): mat3 {
        let cy = Math.cos(yaw);
        let sy = Math.sin(yaw);
        let cp = Math.cos(pitch);
        let sp = Math.sin(pitch);
        let cr = Math.cos(roll);
        let sr = Math.sin(roll);
        let m11 = cr * cp;
        let m12 = - sr * cy + cr * sp * sy;
        let m13 = sr * sy + cr * sp * cy;
        let m21 = sr * cp;
        let m22 = cr * cy + sr * sp * sy;
        let m23 = -cr * sy + sr * sp * cy;
        let m31 = -sp;
        let m32 = cp * sy;
        let m33 = cp * cy;
        return new mat3(
            m11, m12, m13,
            m21, m22, m23,
            m31, m32, m33);
    }
    set(i: number, j: number, value: number): void {
        this.data[i][j] = value;
    }
    get(i: number, j: number): number {
        return this.data[i][j];
    }
    toQuat(): quat {
        throw new Error("Not implemented");
    }
    toMat4(): mat4 {
        return new mat4(
            this.data[0][0], this.data[0][1], this.data[0][2], 0.,
            this.data[1][0], this.data[1][1], this.data[1][2], 0.,
            this.data[2][0], this.data[2][1], this.data[2][2], 0.,
            0., 0., 0., 1.
        );
    }
    static scale(a: mat3, l: number): mat3 {
        return a.scale(l);
    }
    scale(l: number, out?: mat3): mat3 {
        if (!out) {
            out = this.clone();
        }

        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                out.set(i, j, this.get(i, j) * l);
            }
        }
        return out;
    }
    scaleSelf(l: number): mat3 {
        return this.scale(l, this);
    }
    static mul(a: mat3, b: mat3): mat3 {
        return a.mul(b);
    }
    mul(m: mat3, out?: mat3): mat3 {
        if (!out) {
            out = this.clone();
        }

        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                let result = 0.0;
                result += this.data[i][0] * m.data[0][j];
                result += this.data[i][1] * m.data[1][j];
                result += this.data[i][2] * m.data[2][j];
                out.set(i, j, result);
            }
        }
        return out;
    }
    mulSelf(m:mat3): mat3 {
        return this.clone().mul(m, this);
    }
    determinant(): number {
        return this.get(0, 0) * (this.get(1, 1) * this.get(2, 2) - this.get(1, 2) * this.get(2, 1))
            - this.get(0, 1) * (this.get(1, 0) * this.get(2, 2) - this.get(2, 0) * this.get(1, 2))
            + this.get(0, 2) * (this.get(1, 0) * this.get(2, 1) - this.get(2, 0) * this.get(1, 1));
    }
    static solve(m: mat3, v: vec3): vec3 {
        let d = m.determinant();
        let dx = determinant3x3(
            v.x, m.get(0, 1), m.get(0, 2),
            v.y, m.get(1, 1), m.get(1, 2),
            v.z, m.get(2, 1), m.get(2, 2)
        );
        let dy = determinant3x3(
            m.get(0, 0), v.x, m.get(0, 2),
            m.get(1, 0), v.y, m.get(1, 2),
            m.get(2, 0), v.z, m.get(2, 2)
        );
        let dz = determinant3x3(
            m.get(0, 0), m.get(0, 1), v.x,
            m.get(1, 0), m.get(1, 1), v.y,
            m.get(2, 0), m.get(2, 1), v.z
        );
        return new vec3(dx / d, dy / d, dz / d);
    }
    inverse(): mat3 {
        let d = this.determinant();
        let m11 = determinant2x2(
            this.get(1, 1), this.get(1, 2),
            this.get(2, 1), this.get(2, 2)
        );
        let m12 = determinant2x2(
            this.get(0, 2), this.get(0, 1),
            this.get(2, 2), this.get(2, 1)
        );
        let m13 = determinant2x2(
            this.get(0, 1), this.get(0, 2),
            this.get(1, 1), this.get(1, 2)
        );
        let m21 = determinant2x2(
            this.get(1, 2), this.get(1, 0),
            this.get(2, 2), this.get(2, 0)
        );
        let m22 = determinant2x2(
            this.get(0, 0), this.get(0, 2),
            this.get(2, 0), this.get(2, 2)
        );
        let m23 = determinant2x2(
            this.get(0, 2), this.get(0, 0),
            this.get(1, 2), this.get(1, 0)
        );
        let m31 = determinant2x2(
            this.get(1, 0), this.get(1, 1),
            this.get(2, 0), this.get(2, 1)
        );
        let m32 = determinant2x2(
            this.get(0, 1), this.get(0, 0),
            this.get(2, 1), this.get(2, 0)
        );
        let m33 = determinant2x2(
            this.get(0, 0), this.get(0, 1),
            this.get(1, 0), this.get(1, 1)
        );
        return new mat3(
            m11, m12, m13,
            m21, m22, m23,
            m31, m32, m33
        ).scaleSelf(1.0 / d);
    }
    transpose(): mat3 {
        return new mat3(
            this.get(0, 0), this.get(1, 0), this.get(2, 0),
            this.get(0, 1), this.get(1, 1), this.get(2, 1),
            this.get(0, 2), this.get(1, 2), this.get(2, 2)
        );
    }
    transform(point: vec3): vec3 {
        let result = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
            result[i] += this.get(i, 0) * point.x;
            result[i] += this.get(i, 1) * point.y;
            result[i] += this.get(i, 2) * point.z;
        }
        return new vec3(result[0], result[1], result[2]);
    }
    toString(): string {
        return `[
    [${this.get(0, 0).toFixed(4)}, ${this.get(0, 1).toFixed(4)}, ${this.get(0, 2).toFixed(4)}]
    [${this.get(1, 0).toFixed(4)}, ${this.get(1, 1).toFixed(4)}, ${this.get(1, 2).toFixed(4)}]
    [${this.get(2, 0).toFixed(4)}, ${this.get(2, 1).toFixed(4)}, ${this.get(2, 2).toFixed(4)}]
]`;
    }
    // todo:
    static from2DRotation(angle:number, axis:number = 2): mat3 {
        let ca = Math.cos(angle);
        let sa = Math.sin(angle);
        assert(axis >= 0 && axis < 3, "Wrong axis");
        let result = mat3.empty();
        let f = 2 - axis;
        let s = (f + 1) % 3;
        result.set(f, f, ca);
        result.set(f, s, -sa);
        result.set(s, f, sa);
        result.set(s, s, ca);
        return result;
    }
    static fromTRS(t:vec2, r:number, s:vec2):mat3 {
        throw new Error("Not implemented");
    }
    extractMat2(): mat2 {
        return new mat2(
            this.get(0, 0), this.get(0, 1),
            this.get(1, 0), this.get(1, 1)
        );
    }
}