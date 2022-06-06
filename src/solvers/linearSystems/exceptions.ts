


export class InsufficientRankException extends Error {
    constructor(solver: string) {
        super(`${solver} couldn't solve the system`);
    }
}

export class ConvergenseFailureException extends Error {
    constructor(solver: string) {
        super(`${solver} failed to converge`);
    }
}