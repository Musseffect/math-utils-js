import Matrix from "./denseMatrix";
import { assert, clamp, Tolerance, SmallestTolerance } from "./utils";

export default class Vector {
    data: number[];
    constructor(data: number[]) {
        this.data = data;
    }
    toMatrix(rowMajor: boolean = false): Matrix {
        if (rowMajor)
            return new Matrix(this.data, 1, this.data.length);
        return new Matrix(this.data, this.data.length, 1);
    }
    static outer(a: Vector, b: Vector): Matrix {
        let result: Matrix = Matrix.empty(a.size(), b.size());
        for (let j = 0; j < a.size(); j++) {
            for (let i = 0; i < b.size(); i++) {
                result.set(j, i, a.get(j) * b.get(i));
            }
        }
        return result;
    }
    static near(a: Vector, b: Vector, threshold?: number): boolean {
        assert(a.size() == b.size(), "Vectors should have equal sizes");
        if (!threshold)
            threshold = Tolerance;

        for (let i = 0; i < a.size(); i++) {
            if (Math.abs(a.get(i) - b.get(i)) > threshold)
                return false;
        }
        return true;
    }
    static lerp(a: Vector, b: Vector, t: number): Vector {
        return b.sub(a).scaleSelf(t).addSelf(a);
    }
    copy(): Vector {
        return new Vector(this.data.slice());
    }
    clone(): Vector {
        return new Vector(this.data.slice());
    }
    static negate(v: Vector): Vector {
        let result = [];
        for (let el of v.data)
            result.push(-el);
        return new Vector(result);
    }
    static empty(length: number): Vector {
        let data: number[];
        (data = []).length = length;
        data.fill(0);
        return new Vector(data);
    }
    static generate(size: number, f: (i: number) => number) {
        let result = Vector.empty(size);
        for (let i = 0; i < size; ++i)
            result.set(i, f(i));
        return result;
    }
    static dot(a: Vector, b: Vector) {
        let result = 0;
        for (let i = 0; i < a.size(); i++)
            result += a.data[i] * b.data[i];
        return result;
    }
    static add(a: Vector, b: Vector) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] + b.data[i]);
        return new Vector(result);
    }
    static sub(a: Vector, b: Vector) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] - b.data[i]);
        return new Vector(result);
    }
    static mul(a: Vector, b: Vector) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] * b.data[i]);
        return new Vector(result);
    }
    static div(a: Vector, b: Vector) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] / b.data[i]);
        return new Vector(result);
    }
    static scale(a: Vector, s: number) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] * s);
        return new Vector(result);
    }
    addSelf(b: Vector) {
        for (let i = 0; i < this.data.length; i++)
            this.data[i] += b.data[i];
        return this;
    }
    subSelf(b: Vector) {
        for (let i = 0; i < this.data.length; i++)
            this.data[i] -= b.data[i];
        return this;
    }
    scaleSelf(s: number) {
        for (let i = 0; i < this.data.length; i++)
            this.data[i] *= s;
        return this;
    }
    get(i: number): number {
        return this.data[i];
    }
    set(i: number, value: number): void {
        this.data[i] = value;
    }
    size(): number {
        return this.data.length;
    }
    getSubVector(offset: number, length: number) {
        let resultData = new Array(length);
        for (let i = 0; i < length; i++)
            resultData[i] = this.data[offset + i];
        return new Vector(resultData);
    }
    addSubVector(b: Vector, offset: number): Vector {
        assert(this.size() == b.size() + offset, "Vectors should have matching sizes");
        for (let i = 0; i < b.size(); i++)
            this.data[i + offset] += b.get(i);
        return this;
    }
    subSubVector(b: Vector, offset: number): Vector {
        assert(this.size() == b.size() + offset, "Vectors should have matching sizes");
        for (let i = 0; i < b.size(); i++)
            this.data[i + offset] -= b.get(i);
        return this;
    }
    add(b: Vector, dest?: Vector): Vector {
        if (!dest)
            dest = this;
        assert(this.size() == b.size(), "Vectors should be of equal size");
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] + b.data[i];
        return dest;
    }
    sub(b: Vector, dest?: Vector): Vector {
        if (!dest)
            dest = this;
        assert(this.size() == b.size(), "Vectors should be of equal size");
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] - b.data[i];
        return dest;
    }
    mul(b: Vector, dest?: Vector): Vector {
        if (!dest)
            dest = this;
        assert(this.size() == b.size(), "Vectors should be of equal size");
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] * b.data[i];
        return dest;
    }
    scale(b: number, dest?: Vector): Vector {
        if (!dest)
            dest = this;
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] * b;
        return dest;
    }
    l1Norm(): number {
        let result = 0;
        for (let i = 0; i < this.data.length; i++)
            result += Math.abs(this.data[i]);
        return result;
    }
    l2Norm(): number {
        return Math.sqrt(this.squaredLength());
    }
    lInfNorm(): number {
        let result = 0;
        for (let i = 0; i < this.data.length; i++)
            result = Math.max(Math.abs(this.data[i]), result);
        return result;
    }
    lpNorm(p: number): number {
        let result = 0.0;
        for (const value of this.data)
            result += Math.pow(value, p);
        return Math.pow(result, 1 / p);
    }
    squaredLength(): number {
        let result = 0;
        for (let i = 0; i < this.data.length; i++)
            result += this.data[i] * this.data[i];
        return result;
    }
    clamp(min: Vector, max: Vector): void {
        for (let i = 0; i < this.data.length; i++) {
            this.data[i] = clamp(this.data[i], min.data[i], max.data[i]);
        }
    }
    clampScalar(min: number, max: number): void {
        for (let i = 0; i < this.data.length; i++) {
            this.data[i] = clamp(this.data[i], min, max);
        }
    }
    normalize() {
        let l = this.l2Norm();
        if (l > SmallestTolerance) {
            for (let i = 0; i < this.data.length; ++i)
                this.data[i] /= l;
        }
        return this;
    }
    print(fractionDigits: number): string {
        if (!fractionDigits)
            fractionDigits = 4;
        let result = "[ ";
        this.data.forEach((item) => result += item.toFixed(fractionDigits) + " ");
        return result + "]";
    }
    toString(): string {
        let result = "[";
        this.data.forEach((item, i) => { result += `${i != 0 ? ", " : ""}${item}` });
        return result + "]";
    }
}