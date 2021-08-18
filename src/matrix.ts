import { assert, Epsilon } from "./utils";
import vector from "./vector";

export default class matrix {
    w: number;
    h: number;
    data: number[];
    constructor(data: number[], w: number, h: number) {
        this.data = data;
        this.w = w;
        this.h = h;
    }
    static near(a: matrix, b: matrix, threshold?: number): boolean {
        assert(a.w == b.w && a.h == b.h, "Matrices should have equal sizes");
        if (!threshold)
            threshold = Epsilon;

        for (let i = 0; i < a.w; i++) {
            for (let j = 0; j < a.h; j++) {
                if (Math.abs(a.get(j, i) - b.get(j, i)) > threshold)
                    return false;
            }
        }
        return true;
    }
    static random(w: number, h: number): matrix {
        let data = [];
        for (let i = 0; i < w * h; i++)
            data.push(Math.random());
        return new matrix(data, w, h);
    }
    static empty(w: number, h: number): matrix {
        let data: number[];
        (data = []).length = w * h;
        data.fill(0);
        return new matrix(data, w, h);
    }
    copy(): matrix {
        return new matrix(this.data.slice(), this.w, this.h);
    }
    clone(): matrix {
        return new matrix(this.data.slice(), this.w, this.h);
    }
    width(): number {
        return this.w;
    }
    height(): number {
        return this.h;
    }
    static identity(size:number): matrix {
        let data: number[];
        (data = []).length = size * size;
        data.fill(0);
        for (let i = 0; i < size; i++) {
            data[i + i * size] = 1;
        }
        return new matrix(data, size, size);
    }
    static mult(a: matrix, b: matrix): matrix {
        assert(a.w == b.h, "Matrix dimensions aren't compatible");
        let result = matrix.empty(b.w, a.h);
        //for each cell in the result
        for (let j = 0; j < a.h; j++) {
            for (let i = 0; i < b.w; i++) {
                let value = 0;
                for (let k = 0; k < a.w; k++) {
                    value += a.get(j, k) * b.get(k, i);
                }
                result.set(j, i, value);
            }
        }
        return result;
    }
    static postMulVec(a: matrix, b: vector): vector {
        assert(a.w == b.data.length, "Width of matrix isn't compatible with vector's length");
        let result = vector.empty(a.h);
        for (let j = 0; j < a.h; j++) {
            let v = 0;
            for (let i = 0; i < a.w; i++) {
                v += a.get(j, i) * b.get(i);
            }
            result.set(j, v);
        }
        return result;
    }
    static preMulVec(a: matrix, b: vector): vector {
        assert(a.w == b.data.length, "Width of matrix isn't compatible with vector's length");
        let result = vector.empty(a.h);
        for (let i = 0; i < a.w; i++) {
            let v = 0;
            for (let j = 0; j < a.h; j++) {
                v += a.get(j, i) * b.get(j);
            }
            result.set(i, v);
        }
        return result;
    }
    get(row: number, column: number): number {
        return this.data[row * this.w + column];
    }
    // TODO: reverse order of arguments
    set(row: number, column: number, value:number): void {
        this.data[row * this.w + column] = value;
    }
    transpose(): matrix {
        let result = [];
        for (let i = 0; i < this.w; i++) {
            for (let j = 0; j < this.h; j++) {
                result.push(this.data[i + j * this.w]);
            }
        }
        return new matrix(result, this.h, this.w);
    }
    static solve(A: matrix, b: vector): vector {
        assert(A.w == b.data.length, "Width of matrix isn't compatible with vector's length");
        assert(A.w == A.h, "Non-square matrix");
        var rang = b.size();
        var x = vector.empty(rang);
        let epsilon = 0.001
        var indexes = new Array(rang);
        for (var i = 0; i < rang; i++) {
            indexes[i] = i;
        }
        for (var l = 0; l < rang; l++) {
            var max = l;
            for (var i = l + 1; i < rang; i++) {
                if (Math.abs(A.get(indexes[i], l)) > Math.abs(A.get(indexes[max], l)))
                    max = i;
            }
            if (max != l) {
                var temp = indexes[l];
                indexes[l] = indexes[max];
                indexes[max] = temp;
            }
            if (Math.abs(A.get(indexes[l], l)) < epsilon) {
                for (var i = 0; i < rang; i++)
                    x.set(i, 0.0);
                return x;
            }
            for (var i = l + 1; i < rang; i++)
                A.set(indexes[l], i, A.get(indexes[l], i) / A.get(indexes[l], l));
            b.set(indexes[l], b.get(indexes[l]) / A.get(indexes[l], l));
            A.set(indexes[l], l, 1);

            for (var i = l + 1; i < rang; i++) {
                for (var k = l + 1; k < rang; k++)
                    A.set(indexes[i], k, A.get(indexes[i], k) - A.get(indexes[i], l) * A.get(indexes[l], k));
                b.set(indexes[i], b.get(indexes[i]) - A.get(indexes[i], l) * b.get(indexes[l]));
                A.set(indexes[i], l, 0);
            }
        }
        x.set(rang - 1, b.get(indexes[rang - 1]));
        for (var i = rang - 2; i > -1; i--) {
            var k = 0.;
            for (var j = i + 1; j < rang; j++) {
                k = k + A.get(indexes[i], j) * x.get(j);
            }
            x.set(i, b.get(indexes[i]) - k);
        }
        return x;
    }
    inverse(): matrix {
        assert(this.w == this.h, "Non-square matrix");
        let result = this.copy();
        for (let i = 0; i < this.w; i++) {
            let v = vector.empty(this.w);
            v.set(i, 1);
            let column = matrix.solve(this.copy(), v);
            for (let j = 0; j < this.h; j++) {
                result.set(j, i, column.get(j));
            }
        }
        return result;
    }
    print(fractionDigits: number): string {
        if (!fractionDigits)
            fractionDigits = 4;
        let result = "";
        for (let j = 0; j < this.h; j++) {
            result += j > 0 ? "\n| " : "| ";
            for (let i = 0; i < this.w; i++)
                result += this.data[i + j * this.w].toFixed(fractionDigits) + " "
            result += "|";
        }
        return result;
    }
    toString(): string {
        let result = "[";
        for (let j = 0; j < this.h; ++j) {
            if (j != 0)
                result += ",";
            result += "\n\t[";
            for (let i = 0; i < this.w; ++i) {
                if (i != 0) {
                    result += ", "
                }
                result += this.get(j, i).toFixed(4);
            }
            result += "]";
        }
        return result + "\n]";
    }
}