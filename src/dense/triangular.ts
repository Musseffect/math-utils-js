// support both triangular and hessenberg

import AbstractDenseMatrix from "./abstractDenseMatrix";

class TriangularMatrix extends AbstractDenseMatrix {
    _shift: number;
    constructor(data: number[], size: number) {
        super(data, size, size);
        throw "Not implemented";
    }
}