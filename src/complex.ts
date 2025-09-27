import vec2 from "./dense/vec2";

export class complex extends vec2 {
    constructor(x: number, y: number) {
        super(x, y);
    }
    static empty(): complex {
        return new complex(0, 0);
    }
    static real(value: number): complex { return new complex(value, 0.0); }
    static im(value: number): complex { return new complex(0.0, value); }
    public real(): number { return this.x; }
    public im(): number { return this.y; }
    public conjugate(): complex {
        return new complex(this.x, - this.y);
    }
    public arg(): number {
        return Math.atan2(this.y, this.x);
    }
    public inverse(): complex {
        let out = this.conjugate();
        let sl = out.squaredLength();
        out.scaleSelf(1.0 / sl);
        return out;
    }
    public toPolar(): complexPolar {
        return new complexPolar(this.length(), Math.atan2(this.y, this.x));
    }
    static mul(a: complex, b: complex): complex {
        let out = complex.empty();
        out.x = a.x * b.x - a.y * b.y;
        out.y = a.y * b.x + a.x * b.y;
        return out;
    }
    static div(a: complex, b: complex): complex {
        let out = complex.empty();
        let bSquaredLength = b.squaredLength();
        out.x = (a.x * b.x + a.y * b.y) / bSquaredLength;
        out.y = (a.y * b.x - a.x * b.y) / bSquaredLength;
        return out;
    }
    static exp(z: complex): complex {
        return new complexPolar(Math.exp(z.x), z.y).toCartesian();
    }
    static log(z: complex): complex {
        return new complex(Math.log(z.length()), z.arg());
    }
    static pow(a: complex, b: complex): complex {
        let theta = a.arg();
        let lnR = Math.log(a.length());
        let r = Math.exp(b.x * lnR - b.y * theta);
        let angle = b.y * lnR + b.x * theta;
        return new complexPolar(r, angle).toCartesian();
    }
}

export class complexPolar {
    private r: number;
    private angle: number;
    constructor(r: number, angle: number) {
        this.r = r;
        this.angle = angle;
    }
    static empty(): complexPolar {
        return new complexPolar(0, 0);
    }
    public toCartesian(): complex {
        return new complex(this.r * Math.cos(this.angle), this.r * Math.sin(this.angle));
    }
    static mul(a: complexPolar, b: complexPolar): complexPolar {
        return new complexPolar(a.r * b.r, a.angle + b.angle);
    }
    static div(a: complexPolar, b: complexPolar): complexPolar {
        return new complexPolar(a.r / b.r, a.angle - b.angle);
    }
    // todo: wrap angle into [0:2PI) region
    // Todo: add methods
    /*
    public pow(power: number | complex | complexPolar): complexPolar {
        if (power instanceof Number) {
            return new complexPolar(Math.pow(this.r, power), power * this.angle);
        }
        if (power instanceof complex) {
            return new complexPolar();
        }
        return new complexPolar();
    }*/
}