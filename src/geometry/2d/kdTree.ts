import { assert, lerp } from "../../utils";
import vec2 from "../../vec2";


enum Coordinate {
    X = 0,
    Y = 1
};
class Node {
    split: number;
    coord: Coordinate;
    left: Node = null;
    right: Node = null;
};

class SearchNode {
    node: Node;
    minSqrDist: number;
    maxSqrDist: number;
    min: vec2;
    max: vec2;
    constructor(node: Node, min: vec2, max: vec2, sqrDistFunc: (p: vec2) => number) {
        this.node = node;
        this.min = min;
        this.max = max;
        this.maxSqrDist = 0.0;
        this.minSqrDist = Number.POSITIVE_INFINITY;
        for (const point of [new vec2(min.x, min.y), new vec2(min.x, max.y), new vec2(max.x, max.y), new vec2(max.x, min.y)]) {
            let sqrDist = sqrDistFunc(point);
            this.minSqrDist = Math.min(this.minSqrDist, sqrDist);
            this.maxSqrDist = Math.max(this.maxSqrDist, sqrDist);
        }
    }
}

class KDTree {
    points: vec2[];
    min: vec2;
    max: vec2;
    root: Node;
    constructor(points: vec2[]) {
        buildTree(points);
    }
    private buildTree(points: vec2[]) {

    }
    private calcMinSqrDistToNode(point: vec2, node: Node): number {

    }
    private calcMaxSqrDistToNode(point: vec2, node: Node): number {

    }
    public closestPoint(point: vec2): vec2 {
        let searchRadiusSqr = Number.POSITIVE_INFINITY;
        const sqrDistFunc = (p: vec2) => {
            return vec2.squaredDistance(point, p);
        };
        let stack: SearchNode[] = [new SearchNode(this.root, this.min, this.max, sqrDistFunc)];
        searchRadiusSqr = stack[0].maxSqrDist;
        const closestPoint: vec2 = null;
        while (stack.length) {
            let currentNode = stack.pop();
            const node = currentNode.node;
            if (currentNode.minSqrDist > searchRadiusSqr) continue;
            if (node.left == node.right) {
                assert(node.left == null, "Incorrect tree");
                new Error("Not implemented");
            } else {
                let centerCoord = lerp(currentNode.min.getCoord(node.coord), currentNode.max.toArray()[node.coord], currentNode.node.split);
                let leftMax = currentNode.max.clone();
                leftMax.setCoord(centerCoord, node.coord);
                let rightMin = currentNode.min.clone();
                rightMin.setCoord(centerCoord, node.coord);
                let leftNode = new SearchNode(node.left, currentNode.min, leftMax, sqrDistFunc);
                let rightNode = new SearchNode(node.right, rightMin, currentNode.max, sqrDistFunc);
                searchRadiusSqr = Math.min(searchRadiusSqr, leftNode.maxSqrDist, rightNode.maxSqrDist);
                if (leftNode.maxSqrDist < rightNode.maxSqrDist) {
                    stack.push(leftNode);
                    stack.push(rightNode);
                } else {
                    stack.push(rightNode);
                    stack.push(leftNode);
                }
            }
        }
    }
    public closestPointToPrimitive<Primitive>(primitive: Primitive) {


    }
}