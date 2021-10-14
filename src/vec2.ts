import { SmallestEpsilon, vec3 } from ".";
import { Epsilon, lerp, SmallEpsilon } from "./utils";
import vector from "./vector";

export default class vec2 {
    x: number;
    y: number;
    constructor(x: number, y: number) {
        this.x = x;
        this.y = y;
    }
    set(v: vec2): void {
        this.x = v.x;
        this.y = v.y;
    }
    static direction(a: number): vec2 {
        let ca = Math.cos(a);
        let sa = Math.sin(a);
        return new vec2(ca, sa);
    }
    leftNormal(): vec2 {
        return new vec2(-this.y, this.x);
    }
    rightNormal(): vec2 {
        return new vec2(this.y, -this.x);
    }
    toArray(): number[] {
        return [this.x, this.y];
    }
    toVector(): vector {
        return new vector(this.toArray());
    }
    static near(a: vec2, b: vec2, threshold?: number): boolean {
        if (!threshold) {
            threshold = Epsilon;
        }
        return vec2.sub(a, b).lInfnorm() <= threshold;
    }
    static lerp(a: vec2, b: vec2, t: number): vec2 {
        return new vec2(
            lerp(a.x, b.x, t),
            lerp(a.y, b.y, t)
        );
    }
    static empty(): vec2 {
        return new vec2(0., 0.);
    }
    static normalize(v: vec2): vec2 {
        return v.clone().normalize();
    }
    normalize(): vec2 {
        let l = this.l2norm();
        if (l < SmallestEpsilon) {
            this.x = 0;
            this.y = 0;
            return this;
        }
        return this.scaleSelf(1 / l);
    }
    negate(): vec2 {
        this.x = -this.x;
        this.y = -this.y;
        return this;
    }
    static negate(v: vec2): vec2 {
        return new vec2(-v.x, -v.y);
    }
    length(): number {
        return Math.sqrt(this.squaredLength());
    }
    squaredLength(): number {
        return this.x * this.x + this.y * this.y;
    }
    static distance(a: vec2, b: vec2): number {
        return a.sub(b).length();
    }
    static squaredDistance(a: vec2, b: vec2): number {
        return a.sub(b).squaredLength();
    }
    l2norm(): number {
        return Math.sqrt(this.squaredLength());
    }
    l1norm(): number {
        return Math.max(Math.abs(this.x), Math.abs(this.y));
    }
    lInfnorm(): number {
        return Math.abs(this.x) + Math.abs(this.y);
    }
    lpnorm(p: number): number {
        return Math.pow(Math.pow(Math.abs(this.x), p) + Math.pow(Math.abs(this.y), p), 1.0 / p);
    }
    mul(v: vec2, out?: vec2): vec2 {
        if (!out)
            out = this.clone();
        out.x = this.x * v.x;
        out.y = this.y * v.y;
        return out;
    }
    div(v: vec2, out?: vec2): vec2 {
        if (!out)
            out = this.clone();
        out.x = this.x / v.x;
        out.y = this.y / v.y;
        return out;
    }
    add(v: vec2, out?: vec2): vec2 {
        if (!out)
            out = this.clone();
        out.x = this.x + v.x;
        out.y = this.y + v.y;
        return out;
    }
    sub(v: vec2, out?: vec2): vec2 {
        if (!out)
            out = this.clone();
        out.x = this.x - v.x;
        out.y = this.y - v.y;
        return out;
    }
    scale(s: number, out?: vec2): vec2 {
        if (!out)
            out = this.clone();
        out.x = this.x * s;
        out.y = this.y * s;
        return out;
    }
    mulSelf(v: vec2): vec2 {
        return this.mul(v, this);
    }
    divSelf(v: vec2): vec2 {
        return this.div(v, this);
    }
    addSelf(v: vec2): vec2 {
        return this.add(v, this);
    }
    subSelf(v: vec2): vec2 {
        return this.sub(v, this);
    }
    scaleSelf(s: number): vec2 {
        return this.scale(s, this);
    }
    static mul(a: vec2, b: vec2): vec2 {
        return a.mul(b);
    }
    static div(a: vec2, b: vec2): vec2 {
        return a.div(b);
    }
    static add(a: vec2, b: vec2): vec2 {
        return a.add(b);
    }
    static sub(a: vec2, b: vec2): vec2 {
        return a.sub(b);
    }
    static scale(a: vec2, s: number): vec2 {
        return a.scale(s);
    }
    static dot(a: vec2, b: vec2): number {
        return a.x * b.x + a.y * b.y;
    }
    dot(v: vec2): number {
        return vec2.dot(this, v);
    }
    static cross(a: vec2, b: vec2): number {
        return a.x * b.y - a.y * b.x;
    }
    clone(): vec2 {
        return new vec2(this.x, this.y);
    }
    toString(): string {
        return `[${this.x.toFixed(4)}, ${this.y.toFixed(4)}]`;
    }
    apply(op: (a: number) => number): vec2 {
        return new vec2(op(this.x), op(this.y));
    }
    static apply(a: vec2, b: vec2, op: (a: number, b: number) => number): vec2 {
        return new vec2(op(a.x, b.x), op(a.y, b.y));
    }
}