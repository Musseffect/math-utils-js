import { assert } from "../../utils";
import vec2 from "../../vec2";

class PointLocator<Impl> {
    points: vec2[];
    constructor(points: vec2[]) {
        this.points = points;
    }
    public closestPoint(point: vec2): vec2 {
        let sqrDist = Number.POSITIVE_INFINITY;
        let closestPoint: vec2 = null;
        for (const p of this.points) {
            let curSqrDist = vec2.squaredDistance(p, point);
            if (curSqrDist < sqrDist) {
                sqrDist = curSqrDist;
                closestPoint = p;
            }
        }
        return closestPoint;
    }
}


}