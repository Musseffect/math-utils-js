import JSGenerator from "./js";



export default abstract class RandomNumberGenerator {
    abstract randomUnit(): number;
    abstract randomInt(): number;
    abstract random(min: number, max: number): number;
}
