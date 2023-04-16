import Matrix from "../../../denseMatrix";
import { assert, SmallEpsilon } from "../../../utils";
import Vector from "../../../vector";
import { ConvergenseFailureException } from "./exceptions";

const SolverName = "'GaussSeidel'";

export default class GaussSeidel {
    static solve(m: Matrix, rhs: Vector, maxIterations: number, tolerance: number = SmallEpsilon, initialGuess?: Vector) {
        assert(m.width() == m.height(), "Matrix isn't square");
        assert(m.width() == rhs.size(), "Dimensions don't match");
        const rank = rhs.size();
        let result: Vector;
        if (initialGuess) {
            assert(rank == initialGuess.size(), "Initial guess doesn't match system rank");
            result = initialGuess.clone();
        } else {
            result = Vector.empty(rank);
        }
        for (let it = 0; it < maxIterations; ++it) {
            let rhsApprox = Vector.empty(rank);
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
        throw new ConvergenseFailureException(SolverName);
    }
}