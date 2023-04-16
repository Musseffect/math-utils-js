
interface Problem1D {
    f: (x: number) => number;
    min: number,
    max: number,
    solutions: { x: number, y: number }[]
}

const TestFunctions1D: Problem1D[] = [
    {
        // convex unimodal
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
        // non-convex unimodal
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
        // multimodal 1
        f: (x: number) => {
            return Math.sin(x) + Math.sin(10 / 3 * x);
        },
        min: -2.7,
        max: 7.5,
        solutions: [{ x: 5.14573529026, y: -1.89959934915 }]
    },
    {
        // multimodal 2
        f: (x: number) => {
            return -(1.4 - 3 * x) * Math.sin(18 * x);
        },
        min: 0,
        max: 1.2,
        solutions: [{ x: 0.966085803827, y: -1.48907253869 }]
    },
    {
        // multimodal 3
        f: (x: number) => {
            return -x * Math.sin(x);
        },
        min: 0,
        max: 10,
        solutions: [{ x: 7.97866571241, y: -7.91672737159 }]
    }
];


describe('Optimization tests', () => {
    test('1d', () => {

    });
    test('Nd', () => {

    });

})

test('Optimization methods test', () => {


});