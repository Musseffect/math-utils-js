import RandomNumberGenerator from "./generator";

const DefaultSeed = 5489;

const W = 32;
const R = 16;
const M1 = 13;
const M2 = 9;

const FACT = 2.32830643653869628906e-10;

function mat0neg(t: number, v: number) {
    return v ^ (v >> t);
}

function mat0pos(t: number, v: number) {
    return v ^ (v << -t);
}

function mat3neg(t: number, v: number) {
    return v << -t;
}

function mat4neg(t: number, b: number, v: number) {
    return v ^ ((v << -t) & b);
}

export default class WELL512a extends RandomNumberGenerator {
    state: number[];
    stateIdx = 0;
    constructor(seed: number = DefaultSeed) {
        super();
        this.reset(seed);
    }
    private V0(): number {
        return this.state[this.stateIdx];
    }
    private VRm1(): number {
        return this.state[(this.stateIdx + 15) & 0xf];
    }
    private VM1(): number {
        return this.state[(this.stateIdx + M1) & 0xf];
    }
    private VM2(): number {
        return this.state[(this.stateIdx + M2) & 0xf];
    }
    reset(seed: number) {
        this.stateIdx = 0
        this.state = Array(R);
        this.state.fill(0);
        this.state[0] = seed & 0xffffffff;
        if (this.state[0] == 0) {
            this.state[0] = 1 << (W - 1);
        }
    }
    randomUnit(): number {
        return this.randomInt() * FACT;
    }
    randomInt(): number {
        const z0 = this.VRm1();
        const z1 = mat0neg(-16, this.V0()) ^ mat0neg(-15, this.VM1());
        const z2 = mat0pos(11, this.VM2());
        this.state[this.stateIdx] = z1 ^ z2;
        this.state[(this.stateIdx + 15) & 0xf] = mat0neg(-2, z0) ^ mat0neg(-18, z1) ^ mat3neg(-29, z2) ^ mat4neg(-5, 0xda442d24, this.state[this.stateIdx]);
        this.stateIdx = (this.stateIdx + 15) & 0xf;
        return this.state[this.stateIdx]
    }
    random(min: number, max: number): number {
        return this.randomUnit() * (max - min) + min;
    }
}