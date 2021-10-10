import vec2 from "./vec2";

export default class complex extends vec2 {
    constructor(x: number, y: number) {
        super(x, y);
    }
    static polar(r: number, theta: number): complex {
        return new complex(r * Math.cos(theta), r * Math.sin(theta));
    }
    static empty(): complex {
        return new complex(0, 0);
    }
    conjugate(): complex {
        return new complex(this.x, - this.y);
    }
    arg(): number {
        return Math.atan2(this.y, this.x);
    }
    inverse(): complex {
        let out = this.conjugate();
        let sl = out.squaredLength();
        out.scaleSelf(1.0 / sl);
        return out;
    }
    static cMul(a: complex, b: complex): complex {
        let out = complex.empty();
        out.x = a.x * b.x - a.y * b.y;
        out.y = a.y * b.x + a.x * b.y;
        return out;
    }
    static cDiv(a: complex, b: complex): complex {
        let out = complex.empty();
        let bSquaredLength = b.squaredLength();
        out.x = (a.x * b.x + a.y * b.y) / bSquaredLength;
        out.y = (a.y * b.x - a.x * b.y) / bSquaredLength;
        return out;
    }
    static cExp(z: complex): complex {
        return complex.polar(Math.exp(z.x), z.y);
    }
    static cLog(z: complex): complex {
        return new complex(Math.log(z.length()), z.arg());
    }
    static cPow(a: complex, b: complex): complex {
        let theta = a.arg();
        let lnR = Math.log(a.length());
        let r = Math.exp(b.x * lnR - b.y * theta);
        let angle = b.y * lnR + b.x * theta;
        return complex.polar(r, angle);
    }
}