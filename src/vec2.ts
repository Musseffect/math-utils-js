import { Epsilon, lerp, SmallEpsilon } from "./utils";
import vector from "./vector";

export default class vec2 {
    x: number;
    y: number;
    constructor(x: number, y: number) {
        this.x = x;
        this.y = y;
    }
    toArray(): number[] {
        return [this.x, this.y];
    }
    toVector(): vector {
        return new vector(this.toArray());
    }
    static near(a:vec2, b:vec2, threshold?:number):boolean{
        if(!threshold){
            threshold = Epsilon;
        }
        return vec2.sub(a, b).lInfnorm() <= threshold;
    }
    static lerp(a:vec2, b:vec2, t: number):vec2{
        return new vec2(
            lerp(a.x, b.x, t),
            lerp(a.y, b.y, t)
        );
    }
    static empty():vec2{
        return new vec2(0., 0.);
    }
    normalize():vec2{
        let l = this.l2norm();
        if (l < SmallEpsilon) {
            this.x = 0;
            this.y = 0;
            return this;
        }
        return this.scaleSelf(1./l);
    }
    normalized():vec2{
        return this.clone().normalize();
    }
    inverse():vec2 {
        return new vec2(-this.x, -this.y);
    }
    length(): number {
        return Math.sqrt(this.squaredLength());
    }
    squaredLength(): number {
        return this.x * this.x + this.y * this.y;
    }
    l2norm():number{
        return Math.sqrt(this.squaredLength());
    }
    l1norm(): number {
        return Math.max(Math.abs(this.x), Math.abs(this.y));
    }
    lInfnorm(): number {
        return Math.abs(this.x) + Math.abs(this.y);
    }
    lpnorm(p:number): number {
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
    static scale(s: number): vec2 {
        return this.scale(s);
    }
    static dot(a: vec2, b: vec2): number {
        return a.x * b.x + a.y * b.y;
    }
    dot(v:vec2):number{
        return vec2.dot(this, v); 
    }
    static cross(a: vec2, b: vec2): number {
        return a.x * b.y - a.y * b.x;
    }
    clone(): vec2 {
        return new vec2(this.x, this.y);
    }
    toString():string{
        return `[${this.x.toFixed(4)}, ${this.y.toFixed(4)}]`;
    }
}