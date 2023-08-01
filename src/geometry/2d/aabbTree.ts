import { assert } from "../../utils";
import vec2 from "../../vec2";

namespace Polyline {
    class AABBTree {
        points: vec2[];
        constructor(points: vec2[]) {
            this.points = points;
        }
        private checkEmptyTree(): void {
            assert(this.points.length == 0, "empty tree");
        }
        public closestPoint(point: vec2): vec2 {
            this.checkEmptyTree();
        }
        public closestPointAndSegment(point: vec2): { closestPoint: vec2, segmentID: number } {
            this.checkEmptyTree();
        }
        public squaredDistance(point: vec2): number {
            this.checkEmptyTree();
        }
    }
}