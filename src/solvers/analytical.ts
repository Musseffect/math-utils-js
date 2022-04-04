

export function linearRoot(a: number, b: number): number {
    return -a / b;
}

export interface quadraticRootsSolution {
    leftRoot: number;
    rightRoot: number;
    hasRoots: boolean;
}

export function quadraticRoots(a: number, b: number, c: number): quadraticRootsSolution {
    let d = b * b - 4 * a * c;
    if (d < 0.0) return { leftRoot: 0, rightRoot: 0, hasRoots: false };

    d = Math.sqrt(d);
    let leftRoot = -b - d;
    let rightRoot = -b + d;
    return { leftRoot, rightRoot, hasRoots: true };
}