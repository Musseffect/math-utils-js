import RandomNumberGenerator from "./generator";
import JSGenerator from "./js";



export function randomNormalDistr(generator: RandomNumberGenerator = new JSGenerator()) {
    const e1 = generator.randomUnit();
    const e2 = generator.randomUnit();
    return Math.sqrt(-2 * Math.log(e1)) * Math.cos(2 * Math.PI * e2);
}