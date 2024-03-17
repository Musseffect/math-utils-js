

import { assert } from "../../utils";
import Vector from "../../vector";
import {
    OptimizationTestFunction,
    AckleyFunc,
    Beale,
    RastriginFunc,
    RosenbrockFunc,
    SphereFunc,
    GoldsteinPriceFunc,
    BoothFunc,
    BukinFunc6,
    MatyasFunc,
    LeviFunc13,
    HimmelblauFunc,
    ThreeHumpCamelFunc,
    EasomFunc,
    CrossInTrayFunc,
    EggholderFunc,
    HolderTableFunc,
    McCormickFunc,
    SchafferFunc2,
    SchafferFunc4,
    StyblinskiTangFunc
} from "./testFunctions";

interface TestFunction1D {
    name: String;
    f: (x: number) => number;
    min: number,
    max: number,
    solutions: { x: number, y: number }[]
}

const TestFunctions1D: TestFunction1D[] = [
    {
        name: "Convex unimodal",
        f: (x: number) => {
            return 5 + Math.pow(x, 2);
        },
        min: -5,
        max: 5,
        solutions: [{
            x: 0,
            y: 5
        }]
    },
    {
        name: "Non-convex unimodal",
        f: (x: number) => {
            return -(x + Math.sin(x)) * Math.exp(-Math.pow(x, 2));
        },
        min: -10,
        max: 10,
        solutions: [{
            x: 0.67957866002,
            y: -0.824239398476
        }]
    },
    {
        name: "Multimodal 1",
        f: (x: number) => {
            return Math.sin(x) + Math.sin(10 / 3 * x);
        },
        min: -2.7,
        max: 7.5,
        solutions: [{ x: 5.14573529026, y: -1.89959934915 }]
    },
    {
        name: "Multimodal 2",
        f: (x: number) => {
            return -(1.4 - 3 * x) * Math.sin(18 * x);
        },
        min: 0,
        max: 1.2,
        solutions: [{ x: 0.966085803827, y: -1.48907253869 }]
    },
    {
        name: "Multimodal 3",
        f: (x: number) => {
            return -x * Math.sin(x);
        },
        min: 0,
        max: 10,
        solutions: [{ x: 7.97866571241, y: -7.91672737159 }]
    }
];

const TestFunctionsND: OptimizationTestFunction[] = [
    RastriginFunc,
    AckleyFunc,
    SphereFunc,
    RosenbrockFunc,
    Beale,
    GoldsteinPriceFunc,
    BoothFunc,
    BukinFunc6,
    MatyasFunc,
    LeviFunc13,
    HimmelblauFunc,
    ThreeHumpCamelFunc,
    EasomFunc,
    CrossInTrayFunc,
    EggholderFunc,
    HolderTableFunc,
    McCormickFunc,
    SchafferFunc2,
    SchafferFunc4,
    StyblinskiTangFunc
];

describe('Validate optimization test cases', () => {
    test.each(TestFunctionsND)('Test function ND $name', (testFunction: OptimizationTestFunction) => {
        if (typeof testFunction.min.p == "number") {
            expect(testFunction.numDimensions).toBe(-1);
            for (let numDimensions = 2; numDimensions <= 4; ++numDimensions) {
                let p = new Vector(new Array(numDimensions).fill(testFunction.min.p));
                let value = testFunction.f(p);
                expect(value).toBeCloseTo(testFunction.min.value);
            }
            expect(typeof testFunction.searchDomain.min == "number").toBeTruthy()
            expect(typeof testFunction.searchDomain.max == "number").toBeTruthy()
            expect(testFunction.min.p).toBeGreaterThanOrEqual(testFunction.searchDomain.min as number);
            expect(testFunction.min.p).toBeLessThanOrEqual(testFunction.searchDomain.max as number);
        } else {
            expect(testFunction.searchDomain.min instanceof Vector).toBeTruthy();
            expect(testFunction.searchDomain.max instanceof Vector).toBeTruthy();
            const min = testFunction.searchDomain.min as Vector;
            const max = testFunction.searchDomain.max as Vector;
            expect(min.size() == testFunction.numDimensions);
            expect(max.size() == testFunction.numDimensions);
            if (testFunction.min.p instanceof Vector) {
                expect(testFunction.min.p.size()).toBe(testFunction.numDimensions);
                expect(testFunction.f(testFunction.min.p)).toBeCloseTo(testFunction.min.value);
                for (let dimension = 0; dimension < testFunction.numDimensions; ++dimension) {
                    expect(testFunction.min.p.get(dimension)).toBeLessThanOrEqual(max.get(dimension));
                    expect(testFunction.min.p.get(dimension)).toBeGreaterThan(min.get(dimension));
                }
            } else if (Array.isArray(testFunction.min.p)) {
                for (let p of testFunction.min.p) {
                    expect(p.size()).toBe(testFunction.numDimensions);
                    expect(testFunction.f(p)).toBeCloseTo(testFunction.min.value);
                    for (let dimension = 0; dimension < testFunction.numDimensions; ++dimension) {
                        expect(p.get(dimension)).toBeLessThanOrEqual(max.get(dimension));
                        expect(p.get(dimension)).toBeGreaterThan(min.get(dimension));
                    }
                }
            }
        }
    });
});

describe.skip('Optimization tests', () => {
    describe('1d', () => {

    });
    describe('Nd', () => {
        test('Newton', () => {

        });
    });

})

test.skip('Optimization methods test', () => {


});