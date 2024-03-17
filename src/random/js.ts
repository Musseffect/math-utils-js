import RandomNumberGenerator from "./generator";


export default class JSGenerator extends RandomNumberGenerator {
    randomUnit() {
        return Math.random();
    }
}