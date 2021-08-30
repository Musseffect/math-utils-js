import mat3 from "./mat3";
import mat4 from "./mat4";
import matrix from "./matrix";
import { determinant2x2, Epsilon } from "./utils";
import vec2 from "./vec2";

export default class mat2 {
    data: number[][];
    constructor(m11: number, m12: number, m21: number, m22: number){
        this.data = [
            [m11, m12],
            [m21, m22]
        ];
    }
    clone(): mat2 {
        return new mat2(this.data[0][0], this.data[0][1], this.data[1][0], this.data[1][1]);
    }
    set(i: number, j: number, value: number): void {
        this.data[i][j] = value;
    }
    get(i: number, j: number): number {
        return this.data[i][j];
    }
    static near(a: mat2, b: mat2, threshold?: number): boolean {
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
    static empty(): mat2 {
        return new mat2(
            0, 0,
            0, 0);
    }
    static identity(): mat2 {
        return new mat2(
            1, 0,
            0, 1);
    }
    determinant(): number {
        return determinant2x2(this.data[0][0], this.data[0][1], this.data[1][0], this.data[1][1]);
    }
    static scale(a: mat2, l: number): mat2 {
        return a.scale(l);
    }
    scale(l: number, out?: mat2): mat2 {
        if (!out) {
            out = this.clone();
        }
        for (let i = 0; i < 2; i++) {
            for (let j = 0; j < 2; j++) {
                out.set(i, j, this.get(i, j) * l);
            }
        }
        return out;
    }
    scaleSelf(l: number): mat2 {
        return this.scale(l, this);
    }
    inverse(): mat2 {
        let d = this.determinant();
        return new mat2(
            this.data[1][1], -this.data[0][1],
            -this.data[1][0], this.data[0][0]
        ).scaleSelf(1.0 / d);
    }
    toString(): string {
        return `[
    [${this.get(0, 0).toFixed(4)}, ${this.get(0, 1).toFixed(4)}]
    [${this.get(1, 0).toFixed(4)}, ${this.get(1, 1).toFixed(4)}]
]`;
    }
    preMulVec(v: vec2): vec2 {
        return new vec2(v.x * this.data[0][0] + v.y* this.data[1][0], v.x * this.data[0][1] + v.y * this.data[1][1]);
    }
    postMulVec(v: vec2): vec2 {
        return new vec2(v.x * this.data[0][0] + v.y * this.data[0][1], v.x * this.data[1][0] + v.y * this.data[1][1]);
    }
    toMatrix(): matrix {
        let data = [];
        for (let j = 0; j < 4; ++j) {
            for (let i = 0; i < 4; ++i) {
                data.push(this.get(i, j));
            }
        }
        return new matrix(data, 4, 4);
    }
    static rotation(angle: number): mat2 {
        let ca = Math.cos(angle);
        let sa = Math.sin(angle);
        return new mat2(ca, -sa, sa, ca);
    }
    toMat3(): mat3 {
        return new mat3(
            this.data[0][0], this.data[0][1], 0.,
            this.data[1][0], this.data[1][1], 0.,
            0., 0., 1.
        );
    }
    transformPoint2D(p:vec2): vec2 {
        return this.postMulVec(p);
    }
}