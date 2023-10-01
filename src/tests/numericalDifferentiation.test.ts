import Matrix from "../denseMatrix";
import { backwardDifference, centralDifference, forwardDifference, secondOrderDifference } from "../numericalDifferentiation";
import { Epsilon, SmallEpsilon } from "../utils";
import Vector from "../vector";

test("Differentiation", () => {
    let testFunc =
    {
        f: (x: Vector) => {
            return x.get(0) * x.get(0) + Math.sin(3.0 * x.get(0) * x.get(1));
        },
        dfdx: (x: Vector): Vector => {
            const cos = Math.cos(3.0 * x.get(0) * x.get(1));
            return new Vector([2.0 * x.get(0) - 3.0 * x.get(1) * cos, -3.0 * x.get(0) * cos]);
        },
        dfdxdy: (x: Vector): Matrix => {
            const sin = Math.sin(3.0 * x.get(0) * x.get(1));
            const cos = Math.cos(3.0 * x.get(0) * x.get(1));
            let result = Matrix.empty(2, 2);
            result.set(0, 0, 2.0 - 9.0 * x.get(1) * x.get(1) * sin);
            result.set(0, 1, -3.0 * cos - 9.0 * x.get(0) * x.get(1) * sin);
            result.set(1, 0, -3.0 * cos - 9.0 * x.get(0) * x.get(1) * sin);
            result.set(1, 1, -9.0 * x.get(0) * x.get(0) * sin);
            return result;
        }
    };
    let points: Vector[] = [
        new Vector([0, 0]), new Vector([-25, 12]), new Vector([-3, -3]), new Vector([0, 5]),
        new Vector([5, 0]), new Vector([14, 14])];
    const step = SmallEpsilon;
    const tolerance = Epsilon
    for (let point of points) {
        const dfdx = testFunc.dfdx(point);
        const dfdxdy = testFunc.dfdxdy(point);
        expect(Vector.sub(dfdx, forwardDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
        expect(Vector.sub(dfdx, backwardDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
        expect(Vector.sub(dfdx, centralDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
        expect(Matrix.sub(dfdxdy, secondOrderDifference(testFunc.f, point, step)).lInfNorm()).toBeLessThan(tolerance);
    }
});