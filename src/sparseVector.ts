import vector from "./vector";
import { assert, SmallEpsilon } from "./utils";

export default class sparseVector {
    length: number;
    nonZeroElements: number[];
    indices: number[];
    tolerance: number;
    constructor(length: number, nonZeroElements: number[], indices: number[], tolerance: number = SmallEpsilon) {
        this.length = length;
        this.nonZeroElements = nonZeroElements;
        this.indices = indices;
        this.tolerance = tolerance;
    }
    clone(): sparseVector {
        return new sparseVector(this.length, this.nonZeroElements.slice(), this.indices.slice(), this.tolerance);
    }
    static add(v1: sparseVector, v2: sparseVector): sparseVector {
        assert(v1.size() == v2.size(), "Vectors should have the same size");
        if (v1.indices.length == 0) return v2.clone();
        if (v2.indices.length == 0) return v1.clone();
        let ID1 = 0;
        let ID2 = 0;
        let result = new sparseVector(v1.size(), [], [], Math.min(v1.tolerance, v2.tolerance));
        while (ID1 != ID2 && ID1 != v2.indices.length) {
            let currentIndex1 = v1.indices[ID1];
            let currentIndex2 = v2.indices[ID2];
            if (currentIndex1 < currentIndex2) {
                result.indices.push(currentIndex1);
                result.nonZeroElements.push(v1.nonZeroElements[ID1]);
                ID1++;
            } else if (currentIndex2 == currentIndex1) {
                result.indices.push(currentIndex2);
                result.nonZeroElements.push(v2.nonZeroElements[ID2]);
                ID2++;
            } else {
                let value = v1.nonZeroElements[currentIndex1];
                if (Math.abs(value) > result.tolerance) {
                    result.indices.push(currentIndex1);
                    result.nonZeroElements.push(value);
                }
                ID1++;
                ID2++;
            }
        }
        return result;
    }
    static sub(v1: sparseVector, v2: sparseVector): sparseVector {

    }
    static mul(v1: sparseVector, v2: sparseVector): sparseVector {

    }
    static div(v1: sparseVector, v2: sparseVector): sparseVector {

    }
    static negate(v: sparseVector): sparseVector {
        let result = v.clone();
        for (let i = 0; i < v.nonZeroElements.length; ++i)
            result.nonZeroElements[i] = -result.nonZeroElements[i];
        return result;
    }
    static normalize(v: sparseVector): sparseVector {
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
    setTolerance(tolerance: number): void {
        this.tolerance = tolerance;
    }
    getTolerance(): number {
        return this.tolerance;
    }
    static fromVector(array: number[], tolerance: number): sparseVector {
        let nonZeroElements: number[] = [];
        let indices: number[] = [];
        for (let i = 0; i < array.length; ++i) {
            let value = array[i];
            if (value > tolerance) {
                nonZeroElements.push(value);
                indices.push(i);
            }
        }
        let sparse = new sparseVector(array.length, nonZeroElements, indices, tolerance);
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
        let r = this.size();
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
        if (value > this.tolerance) {
            this.indices.splice(l,);
        }
    }
    get(index: number) {
        let l = 0;
        let r = this.size();
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
        let r = this.size();
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
}