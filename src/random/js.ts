import RandomNumberGenerator from "./generator";


export default class JSGenerator extends RandomNumberGenerator {
    randomUnit() {
        return Math.random();
    }
    randomInt() {
        return Math.random() * Number.MAX_SAFE_INTEGER;
    }
    random(min: number, max: number): number {
        return min + (max - min) * Math.random();
    }
}