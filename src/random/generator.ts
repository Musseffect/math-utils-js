import JSGenerator from "./js";



export default abstract class RandomNumberGenerator {
    abstract randomUnit(): number;
    randomInt() {
        return Math.random() * Number.MAX_SAFE_INTEGER;
    }
    random(min: number, max: number): number {
        return min + (max - min) * Math.random();
    }
}
