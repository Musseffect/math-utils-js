


export class SingularMatrixException extends Error {
    constructor(solver: string) {
        super(`${solver} couldn't solve the system due to singular matrix`);
    }
}

export class NotPositiveDefiniteMatrixException extends Error {
    constructor(solver: string) {
        super(`${solver} failed on matrix`);
    }
}

export class ConvergenseFailureException extends Error {
    constructor(solver: string) {
        super(`${solver} failed to converge`);
    }
}

export class InsufficientRankException extends Error {
    constructor(solver: string, rank: number) {
        super(`${solver} couldn't solve the system due to singular matrix of rank ${rank}`)
    }
}
