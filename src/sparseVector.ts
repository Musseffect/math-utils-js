import vector from "./vector";
import { assert, SmallestTolerance } from "./utils";

const DefaultTolerance = SmallestTolerance;

export interface SparseVectorElement {
    index: number;
    value: number;
}

export class SparseVector {
    length: number;
    elements: SparseVectorElement[];
    constructor(length: number, elements: SparseVectorElement[]) {
        this.length = length;
        this.elements = elements;
        elements.sort((a: SparseVectorElement, b: SparseVectorElement): number => {
            return a.index - b.index;
        });
    }
    clone(): SparseVector {
        return new SparseVector(this.length, this.elements.slice());
    }
    static near(v1: SparseVector, v2: SparseVector, tolerance: number): boolean {
        // TODO: rewrite
        let div = SparseVector.sub(v1, v2, tolerance);
        return div.elements.length == 0;
    }
    static empty(size: number): SparseVector {
        return new SparseVector(size, []);
    }
    static binaryOp(v1: SparseVector, v2: SparseVector, op: (a: number, b: number) => number, tolerance: number = 0): SparseVector {
        assert(v1.size() == v2.size(), "Vectors should have the same size");
        if (v1.elements.length == 0) return v2.clone();
        if (v2.elements.length == 0) return v1.clone();
        let it1 = 0;
        let it2 = 0;
        let result = new SparseVector(v1.size(), []);
        while (it1 != v1.elements.length || it2 != v2.elements.length) {
            let currentIndex1 = it1 == v1.elements.length ? v1.size() : v1.elements[it1].index;
            let currentIndex2 = it2 == v2.elements.length ? v2.size() : v2.elements[it2].index;
            let left = 0;
            let right = 0;
            let index = 0;
            if (currentIndex1 <= currentIndex2) {
                index = currentIndex1;
                left = v1.elements[it1].value;
                it1++;
            }
            if (currentIndex2 <= currentIndex1) {
                index = currentIndex2;
                right = v2.elements[it2].value;
                it2++;
            }
            let value = op(left, right);
            if (Math.abs(value) > tolerance)
                result.elements.push({ index, value });
        }
        return result;
    }
    static dot(v1: SparseVector, v2: SparseVector): number {
        let result = SparseVector.mul(v1, v2);
        let value = 0.0;
        for (let element of result.elements)
            value += element.value;
        return value;
    }
    static add(v1: SparseVector, v2: SparseVector, tolerance: number = 0): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a + b, tolerance);
    }
    static sub(v1: SparseVector, v2: SparseVector, tolerance: number = 0): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a - b, tolerance);
    }
    static mul(v1: SparseVector, v2: SparseVector, tolerance: number = 0): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a * b, tolerance);
    }
    static div(v1: SparseVector, v2: SparseVector, tolerance: number = 0): SparseVector {
        return this.binaryOp(v1, v2, (a: number, b: number) => a / b, tolerance);
    }
    static negate(v: SparseVector): SparseVector {
        let result = v.clone();
        for (let i = 0; i < v.elements.length; ++i)
            result.elements[i].value = -result.elements[i].value;
        return result;
    }
    static normalize(v: SparseVector): SparseVector {
        let length = v.l2Norm();
        let result = v.clone();
        for (let i = 0; i < v.elements.length; ++i)
            result.elements[i].value /= length;
        return result;
    }
    l1Norm(): number {
        let result = 0;
        for (let element of this.elements)
            result += Math.abs(element.value);
        return result;
    }
    l2Norm(): number {
        return Math.sqrt(this.squaredLength());
    }
    lInfNorm(): number {
        let result = 0;
        for (let element of this.elements)
            result = Math.max(Math.abs(element.value), result);
        return result;
    }
    squaredLength(): number {
        let result = 0;
        for (let element of this.elements)
            result += element.value * element.value;
        return result;
    }
    static fromVector(array: number[], tolerance: number = DefaultTolerance): SparseVector {
        let elements: SparseVectorElement[] = [];
        for (let index = 0; index < array.length; ++index) {
            let value = array[index];
            if (Math.abs(value) > tolerance)
                elements.push({ index, value });
        }
        let sparse = new SparseVector(array.length, elements);
        return sparse;
    }
    size(): number {
        return this.length;
    }
    toDense(): vector {
        let dense = vector.empty(this.size());
        for (let element of this.elements)
            dense.set(element.index, element.value);
        return dense;
    }
    set(index: number, value: number) {
        let l = 0;
        let r = this.elements.length;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.elements[middle].index < index)
                l = middle + 1;
            else if (this.elements[middle].index > index)
                r = middle;
            else if (value != 0) {
                this.elements[middle].value = value;
                return;
            } else {
                this.elements.splice(middle, 1);
                return;
            }
        }
        if (Math.abs(value) != 0)
            this.elements.splice(l, 0, { value, index });
    }
    get(index: number) {
        let l = 0;
        let r = this.elements.length;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.elements[middle].index < index)
                l = middle + 1;
            else if (this.elements[middle].index > index)
                r = middle;
            else
                return this.elements[middle].value;
        }
        return 0.0;
    }
    isIndexPresent(index: number): boolean {
        let l = 0;
        let r = this.elements.length;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.elements[middle].index < index)
                l = middle + 1;
            else if (this.elements[middle].index > index)
                r = middle;
            else
                return true;
        }
        return false;
    }
    toString(): string {
        let result = `sparse(${this.length})[`;
        for (let i = 0; i < this.elements.length; ++i) {
            result += `${i != 0 ? ", " : ""}${this.elements[i].index}: ${this.elements[i].value}`;
        }
        return result + "]";
    }
}