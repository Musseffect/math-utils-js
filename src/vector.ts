import { assert, clamp } from "./utils";

export default class vector {
    data;
    constructor(data) {
        this.data = data;
    }
    copy(): vector {
        return new vector(this.data.slice());
    }
    clone(): vector {
        return new vector(this.data.slice());
    }
    static empty(length: number): vector {
        let data;
        (data = []).length = length;
        data.fill(0);
        return new vector(data);
    }
    static dot(a, b) {
        let result = 0;
        for (let i = 0; i < a.length(); i++) {
            result += a.data[i] * b.data[i];
        }
        return result;
    }
    static add(a, b) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] + b.data[i]);
        return new vector(result);
    }
    static sub(a, b) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] - b.data[i]);
        return new vector(result);
    }
    static mul(a, b) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] * b.data[i]);
        return new vector(result);
    }
    static scale(a, b) {
        let result = [];
        for (let i = 0; i < a.data.length; i++)
            result.push(a.data[i] * b);
        return new vector(result);
    }
    addSelf(b) {
        for (let i = 0; i < this.data.length; i++)
            this.data[i] += b.data[i];
        return this;
    }
    subSelf(b) {
        for (let i = 0; i < this.data.length; i++)
            this.data[i] -= b.data[i];
        return this;
    }
    scaleSelf(b) {
        for (let i = 0; i < this.data.length; i++)
            this.data[i] *= b;
        return this;
    }
    get(i:number):number {
        return this.data[i];
    }
    set(i:number, value:number):void {
        this.data[i] = value;
    }
    size(): number {
        return this.data.length;
    }
    getSubVector(offset, length) {
        let resultData = new Array(length);
        for (let i = 0; i < length; i++)
            resultData[i] = this.data[offset + i];
        return new vector(resultData);
    }
    addSubVector(b:vector, offset:number):vector {
        assert(this.size() == b.size() + offset, "Vectors should have matching sizes");
        for (let i = 0; i < b.size(); i++)
            this.data[i + offset] += b.get(i);
        return this;
    }
    subSubVector(b:vector, offset:number):vector {
        assert(this.size() == b.size() + offset, "Vectors should have matching sizes");
        for (let i = 0; i < b.size(); i++)
            this.data[i + offset] -= b.get(i);
        return this;
    }
    add(b: vector, dest?: vector): vector {
        if (!dest)
            dest = this;
            assert(this.size() == b.size(), "Vectors should be of equal size");
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] + b.data[i];
        return dest;
    }
    sub(b: vector, dest?: vector): vector {
        if (!dest)
            dest = this;
            assert(this.size() == b.size(), "Vectors should be of equal size");
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] - b.data[i];
        return dest;
    }
    mul(b: vector, dest?: vector): vector {
        if (!dest)
            dest = this;
        assert(this.size() == b.size(), "Vectors should be of equal size");
        for (let i = 0; i < this.data.length; i++)
            dest.data[i] = this.data[i] * b.data[i];
        return dest;
    }
    scale(b: number, dest?: vector): vector {
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
        let result = 0;
        for (let i = 0; i < this.data.length; i++)
            result += this.data[i] * this.data[i];
        return Math.sqrt(result);
    }
    lInfNorm(): number {
        let result = 0;
        for (let i = 0; i < this.data.length; i++)
            result = Math.max(Math.abs(this.data[i]));
        return result;
    }
    norm2(): number {
        let result = 0;
        for (let i = 0; i < this.data.length; i++)
            result += this.data[i] * this.data[i];
        return Math.sqrt(result);
    }
    clamp(min: vector, max: vector): void {
        for (let i = 0; i < this.data.length; i++) {
            this.data[i] = clamp(this.data[i], min.data[i], max.data[i]);
        }
    }
    clampScalar(min: number, max: number): void {
        for (let i = 0; i < this.data.length; i++) {
            this.data[i] = clamp(this.data[i], min, max);
        }
    }
    print(fractionDigits: number): string {
        if (!fractionDigits)
            fractionDigits = 4;
        let result = "[ ";
        this.data.forEach((item) => result += item.toFixed(fractionDigits) + " ");
        return result + "]";
    }
}