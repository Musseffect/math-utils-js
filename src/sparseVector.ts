import vector from "./vector";
import { assert, SmallestTolerance } from "./utils";

const DefaultTolerance = SmallestTolerance;

export default class SparseVector {
    length: number;
    nonZeroElements: number[];
    indices: number[];
    tolerance: number;
    constructor(length: number, nonZeroElements: number[], indices: number[], tolerance: number = DefaultTolerance) {
        this.length = length;
        this.nonZeroElements = nonZeroElements;
        this.indices = indices;
        this.tolerance = tolerance;
    }
    clone(): SparseVector {
        return new SparseVector(this.length, this.nonZeroElements.slice(), this.indices.slice(), this.tolerance);
    }
    static empty(size: number, tolerance: number = DefaultTolerance): SparseVector {
        return new SparseVector(size, [], [], tolerance);
    }
    static binaryOp(v1: SparseVector, v2: SparseVector, op: (a: number, b: number) => number): SparseVector {

        assert(v1.size() == v2.size(), "Vectors should have the same size");
        if (v1.indices.length == 0) return v2.clone();
        if (v2.indices.length == 0) return v1.clone();
        let ID1 = 0;
        let ID2 = 0;
        let result = new SparseVector(v1.size(), [], [], Math.min(v1.tolerance, v2.tolerance));
        while (ID1 != v1.indices.length || ID2 != v2.indices.length) {
            let currentIndex1 = ID1 == v1.indices.length ? v1.size() : v1.indices[ID1];
            let currentIndex2 = ID2 == v2.indices.length ? v2.size() : v2.indices[ID2];
            let left = 0.0;
            let right = 0.0;
            let index = currentIndex1;
            if (currentIndex1 < currentIndex2) {
                left = v1.nonZeroElements[ID1];
                ID1++;
            } else if (currentIndex1 > currentIndex2) {
                index = currentIndex2;
                right = v2.nonZeroElements[ID2];
                ID2++;
            } else {
                left = v1.nonZeroElements[ID1];
                right = v2.nonZeroElements[ID2];
                ID1++;
                ID2++;
            }
            let value = op(left, right);
            if (Math.abs(value) > result.tolerance) {
                result.indices.push(index);
                result.nonZeroElements.push(value);
            }
        }
        return result;
    }
    static dot(v1: SparseVector, v2: SparseVector): number {
        let result = SparseVector.mul(v1, v2);
        let value = 0.0;
        for (let element of result.nonZeroElements)
            value += element;
        return value;
    }
    static add(v1: SparseVector, v2: SparseVector): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a + b);
    }
    static sub(v1: SparseVector, v2: SparseVector): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a - b);
    }
    static mul(v1: SparseVector, v2: SparseVector): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a * b);
    }
    static div(v1: SparseVector, v2: SparseVector): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a / b);
    }
    static negate(v: SparseVector): SparseVector {
        let result = v.clone();
        for (let i = 0; i < v.nonZeroElements.length; ++i)
            result.nonZeroElements[i] = -result.nonZeroElements[i];
        return result;
    }
    static normalize(v: SparseVector): SparseVector {
        let length = v.l2Norm();
        let result = v.clone();
        for (let i = 0; i < v.nonZeroElements.length; ++i)
            result.nonZeroElements[i] = result.nonZeroElements[i] / length;
        return result;
    }
    l1Norm(): number {
        let result = 0;
        for (let element of this.nonZeroElements)
            result += Math.abs(element);
        return result;
    }
    l2Norm(): number {
        return Math.sqrt(this.squaredLength());
    }
    lInfNorm(): number {
        let result = 0;
        for (let element of this.nonZeroElements)
            result = Math.max(Math.abs(element), result);
        return result;
    }
    squaredLength(): number {
        let result = 0;
        for (let element of this.nonZeroElements)
            result += element * element;
        return result;
    }
    setTolerance(tolerance: number = DefaultTolerance): void {
        this.tolerance = tolerance;
    }
    getTolerance(): number {
        return this.tolerance;
    }
    static fromVector(array: number[], tolerance: number = DefaultTolerance): SparseVector {
        let nonZeroElements: number[] = [];
        let indices: number[] = [];
        for (let i = 0; i < array.length; ++i) {
            let value = array[i];
            if (Math.abs(value) > tolerance) {
                nonZeroElements.push(value);
                indices.push(i);
            }
        }
        let sparse = new SparseVector(array.length, nonZeroElements, indices, tolerance);
        return sparse;
    }
    size(): number {
        return this.length;
    }
    toDense(): vector {
        let dense = vector.empty(this.size());
        for (let i = 0; i < this.indices.length; ++i)
            dense.set(this.indices[i], this.nonZeroElements[i]);
        return dense;
    }
    set(index: number, value: number) {
        let l = 0;
        let r = this.indices.length;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.indices[middle] < index)
                l = middle + 1;
            else if (this.indices[middle] > index)
                r = middle;
            else if (Math.abs(value) > this.tolerance) {
                this.nonZeroElements[middle] = value;
                return;
            } else {
                this.nonZeroElements.splice(middle, 1);
                this.indices.splice(middle, 1);
                return;
            }
        }
        if (Math.abs(value) > this.tolerance) {
            this.indices.splice(l, 0, index);
            this.nonZeroElements.splice(l, 0, value);
        }
    }
    get(index: number) {
        let l = 0;
        let r = this.indices.length;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.indices[middle] < index)
                l = middle + 1;
            else if (this.indices[middle] > index)
                r = middle;
            else
                return this.nonZeroElements[middle];
        }
        return 0.0;
    }
    isNonZero(index: number): boolean {
        let l = 0;
        let r = this.indices.length;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.indices[middle] < index)
                l = middle + 1;
            else if (this.indices[middle] > index)
                r = middle;
            else
                return true;
        }
        return false;
    }
    toString(): string {
        let result = `sparse(${this.length})[`;
        for (let i = 0; i < this.indices.length; ++i) {
            result += `${i != 0 ? ", " : ""}${this.indices[i]}: ${this.nonZeroElements[i]}`;
        }
        return result + "]";
    }
}