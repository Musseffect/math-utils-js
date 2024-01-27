import RandomNumberGenerator from "./generator";

// 2^113 period
export default class LFSR113 extends RandomNumberGenerator {
    private x1:number; // should be larger than 1
    private x2:number; // larger than 7
    private x3:number; // larger than 15
    private x4:number; // larger than 127
    constructor(seed: number = 0) {
        super();
        this.reset(seed);
    }
    reset(seed: number) {
        seed = Math.abs(seed);
        this.x1 = Math.max((seed << 1) & 0xFFFFFFFF, 1 << 1);
        this.x2 = Math.max((seed << 3) & 0xFFFFFFFF, 1 << 3);
        this.x3 = Math.max((seed << 4) & 0xFFFFFFFF, 1 << 4);
        this.x4 = Math.max((seed << 5) & 0xFFFFFFFF, 1 << 5);
    }
    randomUnit(): number {
        const m1 = 4294967294;
        const m2 = 4294967288;
        const m3 = 4294967280;
        const m4 = 4294967168;
        const ls1 = 6;
        const ls2 = 2;
        const ls3 = 13;
        const ls4 = 3;
        const rs1 = 13;
        const rs2 = 27;
        const rs3 = 21;
        const rs4 = 12
        let t = (((this.x1 << ls1) ^ this.x1) >> rs1) & 0xFFFFFFFF;
        this.x1 = (((this.x1 & m1) << 18) - t) & 0xFFFFFFFF;
        t = (((this.x2 << ls2) ^ this.x2) >> rs2) & 0xFFFFFFFF;
        this.x2 = (((this.x2 & m2) << 2) - t) & 0xFFFFFFFF;
        t = (((this.x3 << ls3) ^ this.x3) >> rs3) & 0xFFFFFFFF;
        this.x3 = (((this.x3 & m3) << 7) - t) & 0xFFFFFFFF;
        t = (((this.x4 << ls4) ^ this.x4) >> rs4) & 0xFFFFFFFF;
        this.x4 = (((this.x4 & m4) << 13) - t) & 0xFFFFFFFF;
        return ((this.x1 ^ this.x2 ^ this.x3 ^ this.x4) & 0xFFFFFFFF) * 2.3283064365387e-10;
    }
    randomInt(): number {
        throw new Error("Method not implemented.");
    }
    random(min: number, max: number): number {
        throw new Error("Method not implemented.");
    }
}