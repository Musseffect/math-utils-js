import Vector from "../../vector";
import { LineSearch } from "./lineSearch";


export default class WolfeLineSearch extends LineSearch {

    public override step(x: Vector, direction: Vector, initialStep: number = 1.0): number {
        throw new Error("Not implemented");
    }
}