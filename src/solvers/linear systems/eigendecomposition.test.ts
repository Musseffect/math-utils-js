import Matrix from "../../denseMatrix";
import { complex } from "../../complex";

interface TestCase {
    matrix: Matrix;
    eigenvalues: complex[];
    eigendecomposition: {
        Q: Matrix;
        L: Matrix;
    }
}

let testCases: {
    general: TestCase[];
    posDef: TestCase[];
};

(function () {



})();


describe.skip("Eigenvalues extraction", () => {


});