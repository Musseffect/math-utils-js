import matrix from "../../denseMatrix";
import { assert, SmallEpsilon } from "../../utils";
import vector from "../../vector";


// https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
export default class gaussSeidel {
    static solve(m: matrix, rhs: vector, maxIterations: number, tolerance: number = SmallEpsilon, initialGuess?: vector) {
        assert(m.width() == m.height(), "Matrix isn't square");
        assert(m.width() == rhs.size(), "Dimensions don't match");
        const rank = rhs.size();
        let result: vector;
        if (initialGuess) {
            assert(rank == initialGuess.size(), "Initial guess doesn't match system rank");
            result = initialGuess.clone();
        } else {
            result = vector.empty(rank);
        }
        for (let it = 0; it < maxIterations; ++it) {
            let rhsApprox = vector.empty(rank);
            for (let i = 0; i < rank; ++i) {
                let sum = 0.0;
                for (let j = 0; j < i; ++j)
                    sum += m.get(i, j) * result.get(j);
                for (let j = i + 1; j < rank; ++j)
                    sum += m.get(i, j) * result.get(j);
                result.set(i, (rhs.get(i) - sum) / m.get(i, i));
                for (let j = 0; j < rank; ++j)
                    rhsApprox.set(j, rhsApprox.get(j) + m.get(j, i) * result.get(i));
            }
            if (rhsApprox.subSelf(rhs).lInfNorm() < tolerance)
                return result;
        }
        throw new Error("Didn't converge");
    }
}