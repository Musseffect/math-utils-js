import mat3 from "./mat3";
import mat4 from "./mat4";
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

}