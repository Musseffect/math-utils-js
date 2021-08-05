import { Epsilon } from "./utils";

export default class vec2 {
    x: number;
    y: number;
    constructor(x: number, y: number) {
        this.x = x;
        this.y = y;
    }
    static near(a:vec2, b:vec2, threshold?:number):boolean{
        if(!threshold){
            threshold = Epsilon;
        }
        return vec2.sub(a, b).lInfnorm() <= threshold;
    }
    length(): number {
        return Math.sqrt(this.squaredLength());
    }
    squaredLength(): number {
        return this.x * this.x + this.y * this.y;
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
        return this.mul(v, this);
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
    static cross(a: vec2, b: vec2): number {
        return a.x * b.y - a.y * b.x;
    }
    clone(): vec2 {
        return new vec2(this.x, this.y);
    }
}