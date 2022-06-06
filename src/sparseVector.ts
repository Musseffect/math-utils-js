import vector from "./vector";
import { assert, SmallEpsilon } from "./utils";

class sparseVector {
    size: number;
    nonZeroElements: number[];
    indices: number[];
    constructor(size: number, nonZeroElements: number[], indices: number[]) {
        this.size = size;
        this.nonZeroElements = nonZeroElements;
        this.indices = indices;
    }
    static fromVector(array: number[], tolerance: number = SmallEpsilon): sparseVector {
        let nonZeroElements: number[] = [];
        let indices: number[] = [];
        for (let i = 0; i < array.length; ++i) {
            let value = array[i];
            if (value > tolerance) {
                nonZeroElements.push(value);
                indices.push(i);
            }
        }
        return new sparseVector(array.length, nonZeroElements, indices);
    }
    set(index: number, value: number) {
        let l = 0;
        let r = this.size;
        while (l != r) {
            let middle = Math.floor((r + l) / 2);
            if (this.indices[middle] < index)
                l = middle + 1;
            else if (this.indices[middle] > index)
                r = middle;
            else if (value > tolerance) {
                this.nonZeroElements[middle] = value;
                return;
            } else {
                this.nonZeroElements.splice(middle, 1);
                this.indices.splice(middle, 1);
                return;
            }
        }
        if (value > tolerance) {
            this.indices.splice(l,);
        }
    }
    get(index: number) {
        let l = 0;
        let r = this.size;
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
        let r = this.size;
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